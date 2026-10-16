"""
Observation log: everything about the real observations that the stat-binary
files do not record (star identity, magnitudes, coordinates, good-time
windows, faulty mirrors, eLog notes, where the DC / CC reports live).

The log itself is data/observations/observation_log.yaml; this module loads
it into small dataclasses, gives access to each observation's light curves,
renders it as pandas tables and lists what is still missing (``validate``).

Paths in the YAML are resolved against:

- data files (``file``, ``dc_report``, ``conditions_dir``): the observations
  directory, i.e. the YAML's own directory, or ``$KBO_OBSERVATIONS_DIR`` if set;
- ``cc_reports``: ``paths.cc_reports_root`` from the YAML, or
  ``$KBO_CC_REPORTS_ROOT`` if set (the raw CC reports are large and live
  outside the repo).

Absolute paths and ``~`` are accepted everywhere.
"""

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import yaml

from kbo_occultation import PACKAGE_DATA
from kbo_occultation import config

DEFAULT_LOG_PATH = os.path.join(PACKAGE_DATA, "observations", "observation_log.yaml")
OBSERVATIONS_DIR = os.environ.get("KBO_OBSERVATIONS_DIR") or os.path.dirname(DEFAULT_LOG_PATH)


class ObservationLogError(KeyError):
    """A night or observation asked for is not in the log."""

    def __str__(self):
        return self.args[0] if self.args else ""


def _resolve(path, base_dir):
    if path is None:
        return None
    path = os.path.expanduser(path)
    return path if os.path.isabs(path) or base_dir is None else os.path.join(base_dir, path)


@dataclass
class Observation:
    name: str
    night: str
    path: str
    type: str
    ra_deg: Optional[float]
    dec_deg: Optional[float]
    mag_B: Optional[float]
    mag_V: Optional[float]
    spectral_type: Optional[str]
    good_time: Optional[List[float]]
    notes: Optional[str] = None
    manual: Dict = field(default_factory=dict)
    # {telescope: number of panels not correctly focused}: the night's value,
    # updated with this observation's own override.
    faulty_mirrors: Dict = field(default_factory=dict)
    # Channels recording star light (the night's ``channels``).
    channels: List[str] = field(default_factory=list)

    @property
    def is_star(self):
        return self.type == "star"

    @property
    def raw_path(self):
        """The raw stat-binary file the .npz cache is built from."""
        return os.path.splitext(self.path)[0] + ".bin"

    @property
    def has_cache(self):
        return os.path.exists(self.path)

    def time_span(self):
        """
        (t_start, t_stop) of the run in unix seconds, read from the .npz
        cache's t0/dt/n without loading the light curve.
        """
        with np.load(self.path, allow_pickle=False) as npz:
            t0, dt, n = float(npz["t0"]), float(npz["dt"]), int(npz["n"])
        return t0, t0 + (n - 1) * dt

    def magnitude(self, band):
        return {"B": self.mag_B, "V": self.mag_V}[band]

    def lightcurve(self, channel, *, cleaned=True, good_time=False):
        """
        One channel's light curve from the .npz cache (outlier-cleaned by
        default), optionally restricted to the log's good-time window.
        """
        from .photometry import LightCurve, good_time_lightcurve

        lc = LightCurve.from_preprocessed(self.path, channel, cleaned=cleaned)
        return good_time_lightcurve(lc, self.good_time) if good_time else lc

    def lightcurves(self, *, cleaned=True, good_time=False, channels=None):
        """{channel: LightCurve} for the night's channels (or ``channels``)."""
        return {ch: self.lightcurve(ch, cleaned=cleaned, good_time=good_time)
                for ch in (channels or self.channels)}


@dataclass
class Night:
    name: str
    date: str
    dc_report: Optional[str]
    cc_reports: Dict[str, str]
    conditions_dir: Optional[str]
    channels: List[str]
    optical_filter: Optional[str]
    daq: Dict
    dark_run: Optional[str]
    faulty_mirrors: Dict
    elog_notes: Optional[str]
    manual_conditions: Dict
    observations: Dict[str, Observation]

    @staticmethod
    def telescope(channel) -> config.TelescopeConfig:
        """Telescope recorded on ``channel``, from config.ORM_ARRAY."""
        return config.ORM_ARRAY.by_channel(channel)

    @property
    def telescopes(self):
        """Telescope names recording star light this night, from config.ORM_ARRAY."""
        return [self.telescope(ch).name for ch in self.channels]

    @property
    def stars(self):
        return {k: o for k, o in self.observations.items() if o.is_star}

    @property
    def dark(self):
        return self.observations.get(self.dark_run) if self.dark_run else None

    def time_span(self):
        """(t_start, t_stop) unix s covering every cached run of the night."""
        spans = [o.time_span() for o in self.observations.values() if o.has_cache]
        if not spans:
            raise FileNotFoundError(f"night {self.name}: no .npz caches yet (run `kbo-night prepare`)")
        return min(s[0] for s in spans), max(s[1] for s in spans)


@dataclass
class ObservationLog:
    site: Dict
    nights: Dict[str, Night]
    path: Optional[str] = None
    observations_dir: Optional[str] = None
    cc_reports_root: Optional[str] = None

    def all_observations(self):
        return [obs for night in self.nights.values() for obs in night.observations.values()]

    def get(self, night, observation) -> Observation:
        """The observation ``observation`` of night ``night`` (e.g. "2025-12-16", "HD-17316")."""
        try:
            return self.nights[night].observations[observation]
        except KeyError:
            raise ObservationLogError(f"{night!r}/{observation!r} is not in the observation log "
                           f"(nights: {', '.join(self.nights)})") from None

    def select(self, nights=None) -> List[Night]:
        """Nights by name; None, [] or ["all"] means every night."""
        if not nights or list(nights) == ["all"]:
            return list(self.nights.values())
        missing = [n for n in nights if n not in self.nights]
        if missing:
            raise ObservationLogError(f"not in the observation log: {', '.join(missing)} "
                           f"(nights: {', '.join(self.nights)})")
        return [self.nights[n] for n in nights]

    def observation_for_file(self, path):
        """The Observation whose file has the same basename (any extension) as ``path``."""
        stem = os.path.splitext(os.path.basename(path))[0]
        for obs in self.all_observations():
            if os.path.splitext(os.path.basename(obs.path))[0] == stem:
                return obs
        raise ObservationLogError(f"{path!r} is not in the observation log")


def load_observation_log(path=DEFAULT_LOG_PATH) -> ObservationLog:
    with open(path) as f:
        raw = yaml.safe_load(f)
    base_dir = os.environ.get("KBO_OBSERVATIONS_DIR") or os.path.dirname(os.path.abspath(path))
    paths = raw.get("paths") or {}
    cc_root = os.environ.get("KBO_CC_REPORTS_ROOT") or paths.get("cc_reports_root")
    cc_root = os.path.expanduser(cc_root) if cc_root else None

    nights = {}
    for night_name, n in raw["nights"].items():
        night_faulty = n.get("faulty_mirrors") or {}
        channels = n["channels"]
        observations = {}
        for obs_name, o in (n.get("observations") or {}).items():
            observations[obs_name] = Observation(
                name=obs_name,
                night=night_name,
                path=_resolve(o["file"], base_dir),
                type=o["type"],
                ra_deg=o.get("ra_deg"),
                dec_deg=o.get("dec_deg"),
                mag_B=o.get("mag_B"),
                mag_V=o.get("mag_V"),
                spectral_type=o.get("spectral_type"),
                good_time=o.get("good_time"),
                notes=o.get("notes"),
                manual=o.get("manual") or {},
                faulty_mirrors={**night_faulty, **(o.get("faulty_mirrors") or {})},
                channels=channels,
            )
        nights[night_name] = Night(
            name=night_name,
            date=n["date"],
            dc_report=_resolve(n.get("dc_report"), base_dir),
            cc_reports={tel: _resolve(p, cc_root) for tel, p in (n.get("cc_reports") or {}).items()
                        if p is not None},
            conditions_dir=_resolve(n.get("conditions_dir"), base_dir),
            channels=channels,
            optical_filter=n.get("optical_filter"),
            daq=n.get("daq") or {},
            dark_run=n.get("dark_run"),
            faulty_mirrors=night_faulty,
            elog_notes=n.get("elog_notes"),
            manual_conditions=n.get("manual_conditions") or {},
            observations=observations,
        )
    return ObservationLog(site=raw["site"], nights=nights, path=os.path.abspath(path),
                          observations_dir=base_dir, cc_reports_root=cc_root)


def validate(night: Night) -> List[str]:
    """
    Everything still missing for ``night``, one line each: files that are not
    on disk, and log fields that are null (to fill from SIMBAD, the eLog,
    ``kbo-night prepare`` output, ...).
    """
    problems = []
    for obs in night.observations.values():
        if not obs.has_cache and not os.path.exists(obs.raw_path):
            problems.append(f"{obs.name}: neither {os.path.basename(obs.raw_path)} nor its .npz cache found")
        elif not obs.has_cache:
            problems.append(f"{obs.name}: no .npz cache yet (run `kbo-night prepare {night.name}`)")
        # A null good_time means "use the whole run", which is normal for dark runs.
        fields = ("ra_deg", "dec_deg")
        if obs.is_star:
            fields += ("mag_B", "mag_V", "spectral_type", "good_time")
        missing = [f for f in fields if getattr(obs, f) is None]
        if missing:
            problems.append(f"{obs.name}: null {', '.join(missing)}")
    if night.dc_report is None:
        problems.append("dc_report: not set")
    elif not os.path.exists(night.dc_report):
        problems.append(f"dc_report: {night.dc_report} not found")
    for tel, path in night.cc_reports.items():
        if not os.path.exists(path):
            problems.append(f"cc_reports.{tel}: {path} not found")
    if night.conditions_dir and not os.path.isdir(night.conditions_dir) and night.cc_reports:
        problems.append(f"conditions: not extracted yet (run `kbo-night prepare {night.name}`)")
    if night.dark_run and night.dark_run not in night.observations:
        problems.append(f"dark_run: {night.dark_run!r} is not one of the night's observations")
    problems += [f"faulty_mirrors.{t}: null" for t, v in night.faulty_mirrors.items() if v is None]
    problems += [f"manual_conditions.{k}: null" for k, v in night.manual_conditions.items() if v is None]
    if night.elog_notes is None:
        problems.append("elog_notes: null")
    return problems


def nights_table(log: ObservationLog) -> pd.DataFrame:
    """One row per night: channels/telescopes, filter, dark run, mirrors, notes."""
    rows = []
    for night in log.nights.values():
        rows.append({
            "night": night.name,
            "telescopes": ", ".join(night.telescopes),
            "optical_filter": night.optical_filter or "none",
            "n_stars": len(night.stars),
            "dark_run": night.dark_run,
            "cc_reports": ", ".join(night.cc_reports) or "none",
            "faulty_panels": json.dumps(night.faulty_mirrors),
            "elog_notes": night.elog_notes,
        })
    return pd.DataFrame(rows).set_index("night")


def observations_table(log: ObservationLog) -> pd.DataFrame:
    """One row per observation, with the run's UTC start/stop from its .npz cache."""
    from astropy.time import Time

    rows = []
    for obs in log.all_observations():
        start = stop = None
        if obs.has_cache:
            t0, t1 = obs.time_span()
            start, stop = Time([t0, t1], format="unix", scale="utc").isot
        rows.append({
            "night": obs.night,
            "observation": obs.name,
            "type": obs.type,
            "spectral_type": obs.spectral_type,
            "mag_B": obs.mag_B,
            "mag_V": obs.mag_V,
            "ra_deg": obs.ra_deg,
            "dec_deg": obs.dec_deg,
            "start_utc": start,
            "stop_utc": stop,
            "good_time_s": obs.good_time,
            "notes": obs.notes,
        })
    return pd.DataFrame(rows).set_index(["night", "observation"])
