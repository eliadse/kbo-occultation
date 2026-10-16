"""
Per-night preparation of the real observations (``kbo-night prepare DATE``):
from a night's raw files to an analysis-ready night in observation_log.yaml.

``prepare_night`` writes the night's log entry itself, creating it if needed.
Every field that can be derived is filled automatically:

====================================================  ===========================================
field                                                 source
====================================================  ===========================================
observations, file, type, dark_run, daq               ``Spectrum*.bin`` filenames
dc_report                                             ``DCs/<YYYY_MM_DD>/dc_report.pkl`` if present
cc_reports                                            ``CC_M<n>_<YYYY_MM_DD>_*.rep`` under cc_reports_root
channels                                              the data (``detect_channels``)
star ra/dec, mag_B/V, spectral_type, simbad_id        SIMBAD (``observations.simbad``)
dark-run ra/dec                                       CC-REPORT target RA/Dec during the run
faulty_mirrors                                        AMC unfocused panels (night median, per-run overrides)
====================================================  ===========================================

It also builds the ``.npz`` caches (``photometry.preprocess_observation``)
and the CC conditions CSVs. Values are only *filled*: a field already set in
the log (e.g. corrected by hand) is kept unless ``refresh=True``, and the
manual fields (good_time, notes, elog_notes, manual*) are never touched.
Each filled value is marked with an end-of-line comment naming its source.
"""

import datetime as dt
import glob
import os
import re
import warnings
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from kbo_occultation import config
from .cc_reports import amc_unfocused_panels, load_conditions, read_cc_reports, save_conditions
from .log import DEFAULT_LOG_PATH, OBSERVATIONS_DIR, LogFile, Night, append_key, set_value, validate
from .photometry import preprocess_observation

# Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_HD-4335_10002_20002_20251216T212254.bin
FILENAME_RE = re.compile(
    r"^Spectrum_stats_(?P<rate>\d+)MSa_\w*?_(?P<ohm>\d+)Ohm_(?P<mv>\d+)mV_\d+min_"
    r"(?P<name>.+)_\d+_\d+_(?P<ts>\d{8}T\d{6})$")

# A night runs from noon UTC of its date to noon UTC of the next day.
NIGHT_START_HOUR_UTC = 12

# Channel detection: an unconnected DAQ channel only records the electronic
# pedestal, ~0.4-0.7 variance units, identical in every run (2025-12-16 B/C:
# 0.3 % run-to-run spread). A channel with star/NSB light is above
# PEDESTAL_MAX (>= 2.5 even on dark patches) and changes between runs.
PEDESTAL_MAX = 1.5
MIN_RUN_SPREAD = 0.05

# Warn when SIMBAD and the CC target disagree by more than this (deg).
MAX_TARGET_SEPARATION_DEG = 0.5

AUTO_STAR_FIELDS = ("ra_deg", "dec_deg", "mag_B", "mag_V", "spectral_type", "simbad_id")


def _quiet(*_):
    pass


# ───────────────────────────────────────────────────────────
# Discovery from the filenames
# ───────────────────────────────────────────────────────────

def parse_filename(path: str) -> Optional[Dict]:
    """
    Observation name, UTC start (datetime, from the filename timestamp) and DAQ
    settings of a stat-binary file; None if the name does not follow the
    ``Spectrum_stats_<rate>MSa_..._<ohm>Ohm_<mV>mV_<min>min_<name>_<a>_<b>_<YYYYmmddTHHMMSS>``
    pattern.
    """
    m = FILENAME_RE.match(os.path.splitext(os.path.basename(path))[0])
    if m is None:
        return None
    return {
        "name": m["name"],
        "start_utc": dt.datetime.strptime(m["ts"], "%Y%m%dT%H%M%S"),
        "daq": {"sampling": config.standard_sampling,
                "sample_duration_ns": round(1e3 / int(m["rate"])),
                "range_mV": int(m["mv"]), "impedance_ohm": int(m["ohm"])},
    }


def night_files(date: str, data_dir: str = OBSERVATIONS_DIR) -> List[str]:
    """``Spectrum*.bin`` files in ``data_dir`` started during the night of ``date`` (YYYY-MM-DD)."""
    start = dt.datetime.strptime(date, "%Y-%m-%d").replace(hour=NIGHT_START_HOUR_UTC)
    stop = start + dt.timedelta(days=1)
    out = []
    for path in sorted(glob.glob(os.path.join(data_dir, "Spectrum*.bin"))):
        info = parse_filename(path)
        if info and start <= info["start_utc"] < stop:
            out.append(path)
    return sorted(out, key=lambda p: parse_filename(p)["start_utc"])


def find_dc_report(date: str, data_dir: str = OBSERVATIONS_DIR) -> Optional[str]:
    """``DCs/<YYYY_MM_DD>/dc_report.pkl`` relative to ``data_dir``, if it exists."""
    rel = f"DCs/{date.replace('-', '_')}/dc_report.pkl"
    return rel if os.path.exists(os.path.join(data_dir, rel)) else None


def find_cc_reports(date: str, root: Optional[str]) -> Dict[str, str]:
    """
    ``{"m1": path, "m2": path}`` of the night's CC reports, ``CC_M<n>_<YYYY_MM_DD>_*.rep``
    anywhere under ``root``, as paths relative to ``root``. If a telescope has
    several reports that night, the first (earliest) one is taken, with a warning.
    """
    if not root or not os.path.isdir(root):
        return {}
    found = {}
    pattern = os.path.join(root, "**", f"CC_M[0-9]_{date.replace('-', '_')}_*.rep")
    for path in sorted(glob.glob(pattern, recursive=True)):
        tel = "m" + os.path.basename(path)[4]
        if tel in found:
            warnings.warn(f"several CC reports for {tel} on {date}; using {found[tel]}")
            continue
        found[tel] = os.path.relpath(path, root)
    return found


def _observation_text(path, info, indent="      "):
    is_dark = "dark" in info["name"].lower()
    lines = [f"{info['name']}:",
             f"  file: {os.path.splitext(os.path.basename(path))[0]}.npz",
             f"  type: {'dark' if is_dark else 'star'}",
             "  ra_deg: null", "  dec_deg: null"]
    if not is_dark:
        lines += ["  mag_B: null", "  mag_V: null", "  spectral_type: null", "  simbad_id: null",
                  "  good_time: null  # [start, stop] s from the run start; null = whole run"]
    else:
        lines += ["  good_time: null"]
    lines.append("  notes: null")
    return "\n".join(indent + line for line in lines)


def night_skeleton(date: str, data_dir: str = OBSERVATIONS_DIR) -> str:
    """
    YAML text of a new ``nights:`` entry for the night of ``date``, with what
    the filenames tell (runs, types, dark run, DAQ settings) and null for
    everything ``prepare_night`` fills in afterwards.
    """
    files = night_files(date, data_dir)
    if not files:
        raise FileNotFoundError(f"no Spectrum*.bin files of the night of {date} in {data_dir}")
    infos = [parse_filename(p) for p in files]
    darks = [i["name"] for i in infos if "dark" in i["name"].lower()]
    daq = infos[0]["daq"]
    lines = [
        f'  "{date}":',
        f'    date: "{date}"',
        "    dc_report: null",
        "    cc_reports: null",
        f"    conditions_dir: conditions/{date}",
        "    channels: null",
        "    optical_filter: null",
        "    daq: {" + ", ".join(f"{k}: {v}" for k, v in daq.items()) + "}",
        f"    dark_run: {darks[0] if darks else 'null'}",
        "    faulty_mirrors:",
        "      m1: null",
        "      m2: null",
        "      lst1: null  # not in the MAGIC CC reports",
        "    elog_notes: null  # TODO(user): weather/eLog remarks for the night",
        "    observations:",
    ]
    lines += [_observation_text(p, i) for p, i in zip(files, infos)]
    return "\n".join(lines) + "\n"


def _add_missing_runs(lf: LogFile, date: str, changes: List[str]):
    """Add to an existing night the runs of its date that are on disk but not in the log."""
    entry = lf.nights[date]
    listed = {os.path.splitext(o["file"])[0] for o in (entry.get("observations") or {}).values()}
    for path in night_files(date, lf.base_dir):
        stem = os.path.splitext(os.path.basename(path))[0]
        if stem in listed:
            continue
        info = parse_filename(path)
        new = lf.load_mapping(_observation_text(path, info, indent=""))
        append_key(entry["observations"], info["name"], new[info["name"]])
        changes.append(f"{info['name']}: new run added")


# ───────────────────────────────────────────────────────────
# Derived values
# ───────────────────────────────────────────────────────────

def channel_levels(night: Night) -> pd.DataFrame:
    """Median outlier-cleaned variance per run (rows) and channel A/B/C (columns), from the caches."""
    rows = {}
    for obs in night.observations.values():
        if obs.has_cache:
            with np.load(obs.path, allow_pickle=False) as npz:
                rows[obs.name] = {ch: float(np.median(npz[f"var_{ch}_clean"])) for ch in "ABC"}
    return pd.DataFrame.from_dict(rows, orient="index")


def detect_channels(levels: pd.DataFrame) -> List[str]:
    """
    Channels that recorded light, from ``channel_levels``: median level above
    PEDESTAL_MAX, or run-to-run spread (max/min - 1) above MIN_RUN_SPREAD.
    """
    lit = []
    for ch in levels.columns:
        v = levels[ch].to_numpy()
        spread = v.max() / v.min() - 1.0 if len(v) > 1 and v.min() > 0 else 0.0
        if np.median(v) > PEDESTAL_MAX or spread > MIN_RUN_SPREAD:
            lit.append(ch)
    return lit


def _site_reports(conditions):
    return conditions.get("m2") or next(iter(conditions.values()), {})


def cc_target(conditions, t_start, t_stop):
    """Median CC-REPORT target (RA, Dec) in deg during [t_start, t_stop], or None."""
    cc = _site_reports(conditions).get("cc")
    if cc is None or cc.empty or "target_ra_deg" not in cc:
        return None
    sel = cc[(cc["time"] >= t_start) & (cc["time"] <= t_stop)]
    sel = sel[(sel["target_ra_deg"] != 0) | (sel["target_dec_deg"] != 0)].dropna(subset=["target_ra_deg"])
    if not len(sel):
        return None
    return round(float(sel["target_ra_deg"].median()), 2), round(float(sel["target_dec_deg"].median()), 2)


def _separation_deg(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


# ───────────────────────────────────────────────────────────
# Caches and conditions
# ───────────────────────────────────────────────────────────

def build_caches(night: Night, *, force: bool = False, log: Callable = print) -> List[str]:
    """Write the missing (all, with ``force``) ``.npz`` caches of ``night``. Returns the paths written."""
    written = []
    for obs in night.observations.values():
        if obs.has_cache and not force:
            continue
        if not os.path.exists(obs.raw_path):
            log(f"  {obs.name}: {os.path.basename(obs.raw_path)} not found, no cache built")
            continue
        log(f"  {obs.name}: pre-processing {os.path.basename(obs.raw_path)} ...")
        written.append(preprocess_observation(obs.raw_path, out_path=obs.path))
    return written


def _covers(reports, t_start, t_stop):
    cc = reports.get("cc")
    return cc is not None and not cc.empty and cc["time"].min() <= t_start and cc["time"].max() >= t_stop


def extract_conditions(night: Night, *, force: bool = False, pad_s: float = 1800.0,
                       log: Callable = print) -> Dict[str, Dict]:
    """
    Condense the night's CC ``.rep`` reports into ``<conditions_dir>/<tel>_<report>.csv``,
    restricted to the night's observing span +- ``pad_s``. Existing CSVs are
    reused unless ``force`` or they don't cover every run (e.g. runs added
    later). Returns ``{telescope: {report: DataFrame}}``.
    """
    if not night.cc_reports or night.conditions_dir is None:
        log("  no CC reports found, skipping the conditions")
        return {}
    cached = load_conditions(night.conditions_dir) if os.path.isdir(night.conditions_dir) else {}
    span = night.time_span()
    out = {}
    for tel, rep_path in night.cc_reports.items():
        if tel in cached and not force and _covers(cached[tel], *span):
            out[tel] = cached[tel]
            continue
        if not os.path.exists(rep_path):
            log(f"  {tel}: {rep_path} not found, skipping")
            continue
        log(f"  {tel}: parsing {os.path.basename(rep_path)} ...")
        reports = read_cc_reports(rep_path, span[0] - pad_s, span[1] + pad_s)
        save_conditions(reports, night.conditions_dir, tel)
        log("    " + ", ".join(f"{k}={len(v)}" for k, v in reports.items()))
        out[tel] = reports
    return out


# ───────────────────────────────────────────────────────────
# Filling the log
# ───────────────────────────────────────────────────────────

class _Filler:
    """Fill-only (or, with refresh, overwrite) edits of one night entry, recorded as changes."""

    def __init__(self, refresh: bool):
        self.refresh = refresh
        self.changes: List[str] = []

    def __call__(self, cm, key, value, source, label):
        if value is None:
            return
        old = cm.get(key)
        if old is not None and not self.refresh:
            return
        if _same(old, value):
            return
        if key in cm:
            set_value(cm, key, value, source)
        else:
            append_key(cm, key, value, source)
        self.changes.append(f"{label}{key}: {_fmt(old)} -> {_fmt(value)}  [{source}]")


def _same(a, b):
    if isinstance(a, float) or isinstance(b, float):
        try:
            return abs(float(a) - float(b)) < 1e-9
        except (TypeError, ValueError):
            return False
    if isinstance(b, (list, dict)):
        return a is not None and (list(a) == list(b) if isinstance(b, list) else dict(a) == dict(b))
    return a == b


def _fmt(v):
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(map(str, v)) + "]"
    if hasattr(v, "items"):
        return "{" + ", ".join(f"{k}: {x}" for k, x in v.items()) + "}"
    return str(v)


def _flow_list(lf, items):
    return lf.load_mapping("x: [" + ", ".join(items) + "]")["x"]


def _block_map():
    from ruamel.yaml.comments import CommentedMap
    return CommentedMap()


def _drop_unlit_mirrors(entry, fill: _Filler):
    """Remove null faulty_mirrors entries of telescopes that recorded no light this night."""
    faulty, channels = entry.get("faulty_mirrors"), entry.get("channels")
    if not faulty or not channels:
        return
    lit = {config.ORM_ARRAY.by_channel(ch).dc_key for ch in channels}
    for tel in [t for t, v in faulty.items() if v is None and t not in lit]:
        del faulty[tel]
        fill.changes.append(f"faulty_mirrors.{tel}: removed (no light on its channel)")


def _fill_night(lf: LogFile, date: str, fill: _Filler, *, force: bool, pad_s: float,
                query_stars: Callable, log: Callable):
    entry = lf.nights[date]

    # files on disk
    fill(entry, "dc_report", find_dc_report(date, lf.base_dir), "found", "")
    found = find_cc_reports(date, lf.cc_root)
    if found:
        if entry.get("cc_reports") is None:
            set_value(entry, "cc_reports", _block_map())
        for tel, rel in found.items():
            fill(entry["cc_reports"], tel, rel, "found", "cc_reports.")

    # caches and channels
    night = lf.night(date)
    log(f"── {date}: caches")
    build_caches(night, force=force, log=log)
    levels = channel_levels(night)
    if not levels.empty:
        lit = detect_channels(levels)
        log("  median level per channel (variance units):\n"
            + "\n".join("    " + line for line in levels.round(3).to_string().splitlines()))
        if lit:
            fill(entry, "channels", _flow_list(lf, lit), "data", "")
        else:
            log("  no channel above the pedestal: channels left unset")
    _drop_unlit_mirrors(entry, fill)

    # conditions
    night = lf.night(date)
    log(f"── {date}: conditions")
    conditions = extract_conditions(night, force=force, pad_s=pad_s, log=log)

    # coordinates and magnitudes
    observations = entry["observations"]
    stars = [o.name for o in night.observations.values()
             if o.is_star and (fill.refresh or any(observations[o.name].get(k) is None for k in AUTO_STAR_FIELDS))]
    simbad = query_stars(stars) if stars else {}
    for obs in night.observations.values():
        cm = observations[obs.name]
        target = cc_target(conditions, *obs.time_span()) if obs.has_cache and conditions else None
        if obs.is_star:
            info = simbad.get(obs.name)
            if obs.name in simbad and info is None:
                log(f"  {obs.name}: not resolved by SIMBAD, coordinates/magnitudes left unset")
            for key in AUTO_STAR_FIELDS:
                fill(cm, key, (info or {}).get(key), "SIMBAD", f"{obs.name}.")
            ra, dec = cm.get("ra_deg"), cm.get("dec_deg")
            if target and ra is not None and dec is not None:
                sep = _separation_deg(ra, dec, *target)
                if sep > MAX_TARGET_SEPARATION_DEG:
                    log(f"  WARNING {obs.name}: SIMBAD position is {sep:.2f} deg from the CC target {target}")
        elif target:
            fill(cm, "ra_deg", target[0], "CC target", f"{obs.name}.")
            fill(cm, "dec_deg", target[1], "CC target", f"{obs.name}.")

    # faulty mirrors: night median, per-run override where a run differs
    faulty = entry.get("faulty_mirrors")
    for tel, reports in conditions.items():
        counts = {o.name: amc_unfocused_panels(reports.get("amc"), *o.time_span())
                  for o in night.observations.values() if o.has_cache}
        known = [c for c in counts.values() if c is not None]
        log(f"  {tel} AMC unfocused panels per run: " + ", ".join(f"{k}={v}" for k, v in counts.items()))
        if not known:
            continue
        if faulty is None:
            set_value(entry, "faulty_mirrors", _block_map())
            faulty = entry["faulty_mirrors"]
        if faulty.get(tel) is not None and not fill.refresh:
            continue   # night value set by hand: leave its per-run values alone too
        night_value = int(round(float(np.median(known))))
        fill(faulty, tel, night_value, "AMC", "faulty_mirrors.")
        for name, count in counts.items():
            cm = observations[name]
            per_run = cm.get("faulty_mirrors")
            if count is None or count == night_value:
                if fill.refresh and per_run is not None and tel in per_run:
                    del per_run[tel]
                continue
            if per_run is None:
                append_key(cm, "faulty_mirrors", _block_map())
                per_run = cm["faulty_mirrors"]
            fill(per_run, tel, count, "AMC", f"{name}.faulty_mirrors.")


def prepare_night(date: str, *, log_path: str = DEFAULT_LOG_PATH, force: bool = False,
                  refresh: bool = False, dry_run: bool = False, pad_s: float = 1800.0,
                  query_stars: Optional[Callable] = None, log: Callable = print) -> List[str]:
    """
    Make the night of ``date`` (YYYY-MM-DD) ready for analysis and write its
    entry in the log at ``log_path``: create it from the ``Spectrum*.bin``
    files if it is not there, add runs that are not listed yet, build the
    caches and conditions, and fill every derivable field (see the module
    docstring). ``force`` rebuilds the caches / conditions; ``refresh``
    overwrites already-set automatic fields; ``dry_run`` prints the changes
    without writing the log (caches and conditions are still built).
    ``query_stars`` replaces the SIMBAD lookup (default
    ``observations.simbad.query_stars``). Returns ``validate()`` of the
    result: what is still missing.
    """
    if query_stars is None:
        from .simbad import query_stars
    lf = LogFile(log_path)
    fill = _Filler(refresh)
    if date not in lf.nights:
        lf.append_night(night_skeleton(date, lf.base_dir))
        fill.changes.append("new night")
    else:
        _add_missing_runs(lf, date, fill.changes)

    _fill_night(lf, date, fill, force=force, pad_s=pad_s, query_stars=query_stars, log=log)

    log(f"── {date}: observation log ({lf.path})")
    for c in fill.changes:
        log(f"  {c}")
    if not fill.changes:
        log("  no changes")
    elif dry_run:
        log("  (dry run: not written)")
    else:
        lf.save()

    problems = validate(lf.night(date))
    log(f"── {date}: still missing")
    for p in problems:
        log(f"  {p}")
    if not problems:
        log("  nothing")
    return problems
