"""
Per-night preparation of the real observations, the step between copying a
night's files and analysing them (``kbo-night scaffold`` / ``kbo-night prepare``):

- ``scaffold_night``: the YAML block to paste into observation_log.yaml for a
  new night, built from the night's ``Spectrum*.bin`` filenames;
- ``prepare_night``: build the missing ``.npz`` caches from the raw ``.bin``
  files (``photometry.preprocess_observation``), condense the CC ``.rep``
  reports into conditions CSVs (``extract_conditions``), and report the
  per-run AMC unfocused panels and what is still missing from the log.

Everything is idempotent: existing caches are kept unless ``force=True``.
"""

import datetime as dt
import glob
import os
import re
from typing import Callable, Dict, List, Optional

from kbo_occultation import config
from .cc_reports import amc_unfocused_panels, load_conditions, read_cc_reports, save_conditions
from .log import OBSERVATIONS_DIR, Night, validate
from .photometry import preprocess_observation

# Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_HD-4335_10002_20002_20251216T212254.bin
FILENAME_RE = re.compile(
    r"^Spectrum_stats_(?P<rate>\d+)MSa_\w*?_(?P<ohm>\d+)Ohm_(?P<mv>\d+)mV_\d+min_"
    r"(?P<name>.+)_\d+_\d+_(?P<ts>\d{8}T\d{6})$")

# A night runs from noon UTC of its date to noon UTC of the next day.
NIGHT_START_HOUR_UTC = 12


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


def scaffold_night(date: str, data_dir: str = OBSERVATIONS_DIR) -> str:
    """
    Text of a ``nights:`` entry for observation_log.yaml covering the night of
    ``date``, for you to paste and complete. Filled from the filenames: the
    observation names and files, ``type: dark`` (and ``dark_run``) for names
    containing "dark", the DAQ settings, and ``dc_report`` if
    ``DCs/<YYYY_MM_DD>/dc_report.pkl`` exists. Everything else is
    ``null  # TODO(user)``: coordinates, magnitudes and spectral types from
    SIMBAD, good-time windows, CC report paths, faulty mirrors, eLog notes.
    """
    files = night_files(date, data_dir)
    if not files:
        raise FileNotFoundError(f"no Spectrum*.bin files of the night of {date} in {data_dir}")
    infos = [parse_filename(p) for p in files]
    darks = [i["name"] for i in infos if "dark" in i["name"].lower()]
    dc_rel = f"DCs/{date.replace('-', '_')}/dc_report.pkl"
    dc_line = (f"dc_report: {dc_rel}" if os.path.exists(os.path.join(data_dir, dc_rel))
               else f"dc_report: null  # TODO(user): {dc_rel}")
    daq = infos[0]["daq"]
    todo = "  # TODO(user)"

    lines = [
        f'  "{date}":',
        f'    date: "{date}"',
        f"    {dc_line}",
        "    cc_reports:  # relative to paths.cc_reports_root",
        f"      m1: null{todo}",
        f"      m2: null{todo}",
        f"    conditions_dir: conditions/{date}",
        f"    channels: [A, B, C]  # TODO(user): drop channels without star light",
        "    optical_filter: null",
        "    daq: {" + ", ".join(f"{k}: {v}" for k, v in daq.items()) + "}",
        f"    dark_run: {darks[0] if darks else 'null'}",
        "    faulty_mirrors:  # `kbo-night prepare` prints the AMC values",
        f"      m1: null{todo}",
        f"      m2: null{todo}",
        "      lst1: null  # not in the MAGIC CC reports",
        f"    elog_notes: null{todo}",
        "    observations:",
    ]
    for path, info in zip(files, infos):
        is_dark = "dark" in info["name"].lower()
        lines += [
            f"      {info['name']}:",
            f"        file: {os.path.splitext(os.path.basename(path))[0]}.npz",
            f"        type: {'dark' if is_dark else 'star'}",
            f"        ra_deg: null{todo}",
            f"        dec_deg: null{todo}",
        ]
        if is_dark:
            lines += ["        mag_B: null", "        mag_V: null", "        spectral_type: null",
                      "        good_time: null"]
        else:
            lines += [f"        mag_B: null{todo}", f"        mag_V: null{todo}",
                      f"        spectral_type: null{todo}",
                      f"        good_time: null{todo}: [start, stop] s from the run start"]
        lines.append("        notes: null")
    return "\n".join(lines) + "\n"


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


def extract_conditions(night: Night, *, force: bool = False, pad_s: float = 1800.0,
                       log: Callable = print) -> Dict[str, Dict]:
    """
    Condense the night's CC ``.rep`` reports into ``<conditions_dir>/<tel>_<report>.csv``,
    restricted to the night's observing span +- ``pad_s``. Telescopes whose CSVs
    already exist are only re-read from the CSVs unless ``force``. Returns
    ``{telescope: {report: DataFrame}}``.
    """
    if not night.cc_reports or night.conditions_dir is None:
        log("  no CC reports listed, skipping the conditions")
        return {}
    cached = load_conditions(night.conditions_dir) if os.path.isdir(night.conditions_dir) else {}
    out = {}
    t_start = t_stop = None
    for tel, rep_path in night.cc_reports.items():
        if tel in cached and not force:
            out[tel] = cached[tel]
            continue
        if not os.path.exists(rep_path):
            log(f"  {tel}: {rep_path} not found, skipping")
            continue
        if t_start is None:
            t_start, t_stop = night.time_span()
            t_start, t_stop = t_start - pad_s, t_stop + pad_s
        log(f"  {tel}: parsing {os.path.basename(rep_path)} ...")
        reports = read_cc_reports(rep_path, t_start, t_stop)
        save_conditions(reports, night.conditions_dir, tel)
        log("    " + ", ".join(f"{k}={len(v)}" for k, v in reports.items()))
        out[tel] = reports
    return out


def prepare_night(night: Night, *, force: bool = False, pad_s: float = 1800.0,
                  log: Callable = print) -> List[str]:
    """
    Make ``night`` ready for analysis: build the ``.npz`` caches, extract the
    CC conditions, print the per-run AMC unfocused panels (undefined + error;
    copy them into the log's ``faulty_mirrors``) and return ``validate(night)``,
    the list of what is still missing.
    """
    log(f"── {night.name}: caches")
    build_caches(night, force=force, log=log)
    log(f"── {night.name}: conditions")
    conditions = extract_conditions(night, force=force, pad_s=pad_s, log=log)
    for tel, reports in conditions.items():
        amc = reports.get("amc")
        counts = {obs.name: amc_unfocused_panels(amc, *obs.time_span())
                  for obs in night.observations.values() if obs.has_cache}
        log(f"  {tel} unfocused panels per run (log value: {night.faulty_mirrors.get(tel)}): "
            + ", ".join(f"{k}={v}" for k, v in counts.items()))
    problems = validate(night)
    log(f"── {night.name}: still missing")
    for p in problems:
        log(f"  {p}")
    if not problems:
        log("  nothing")
    return problems
