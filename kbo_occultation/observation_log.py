"""
Observation log: everything about the real observations that the stat-binary
files do not record (star identity, magnitudes, coordinates, good-time
windows, faulty mirrors, eLog notes, where the DC / CC reports live).

The log itself is data/observations/observation_log.yaml; this module only
loads it into small dataclasses and renders it as pandas tables.
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

OBSERVATIONS_DIR = os.path.join(PACKAGE_DATA, "observations")
DEFAULT_LOG_PATH = os.path.join(OBSERVATIONS_DIR, "observation_log.yaml")


def _resolve(path, base_dir):
    if path is None:
        return None
    return path if os.path.isabs(path) else os.path.join(base_dir, path)


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

    @property
    def is_star(self):
        return self.type == "star"

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

    @property
    def telescopes(self):
        """Telescope names recording star light this night, from config.ORM_ARRAY."""
        return [config.ORM_ARRAY.by_channel(ch).name for ch in self.channels]

    @property
    def stars(self):
        return {k: o for k, o in self.observations.items() if o.is_star}

    @property
    def dark(self):
        return self.observations.get(self.dark_run) if self.dark_run else None


@dataclass
class ObservationLog:
    site: Dict
    nights: Dict[str, Night]

    def all_observations(self):
        return [obs for night in self.nights.values() for obs in night.observations.values()]

    def observation_for_file(self, path):
        """The Observation whose file has the same basename (any extension) as ``path``."""
        stem = os.path.splitext(os.path.basename(path))[0]
        for obs in self.all_observations():
            if os.path.splitext(os.path.basename(obs.path))[0] == stem:
                return obs
        raise KeyError(f"{path!r} is not in the observation log")


def load_observation_log(path=DEFAULT_LOG_PATH) -> ObservationLog:
    with open(path) as f:
        raw = yaml.safe_load(f)
    base_dir = os.path.dirname(os.path.abspath(path))

    nights = {}
    for night_name, n in raw["nights"].items():
        observations = {}
        for obs_name, o in n["observations"].items():
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
            )
        nights[night_name] = Night(
            name=night_name,
            date=n["date"],
            dc_report=_resolve(n.get("dc_report"), base_dir),
            cc_reports={tel: _resolve(p, base_dir) for tel, p in (n.get("cc_reports") or {}).items()},
            conditions_dir=_resolve(n.get("conditions_dir"), base_dir),
            channels=n["channels"],
            optical_filter=n.get("optical_filter"),
            daq=n.get("daq") or {},
            dark_run=n.get("dark_run"),
            faulty_mirrors=n.get("faulty_mirrors") or {},
            elog_notes=n.get("elog_notes"),
            manual_conditions=n.get("manual_conditions") or {},
            observations=observations,
        )
    return ObservationLog(site=raw["site"], nights=nights)


def _count_panels(entry):
    if not entry:
        return None
    return len(entry.get("error_panels", [])) + len(entry.get("undefined_panels", []))


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
            "faulty_panels": json.dumps({tel: _count_panels(v) for tel, v in night.faulty_mirrors.items()}),
            "elog_notes": night.elog_notes,
        })
    return pd.DataFrame(rows).set_index("night")


def observations_table(log: ObservationLog) -> pd.DataFrame:
    """One row per observation, with the run's UTC start/stop from its .npz cache."""
    from astropy.time import Time

    rows = []
    for obs in log.all_observations():
        start = stop = None
        if os.path.exists(obs.path):
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
