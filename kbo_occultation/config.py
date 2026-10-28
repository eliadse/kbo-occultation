# kbo_occultation/config.py

from dataclasses import dataclass
from typing import Tuple

import numpy as np

# ── Noise-removal constants (single source of truth) ──────────────────────
# Butterworth high-pass cutoff used by the "highpass" correction everywhere.
HIGHPASS_CUTOFF_HZ = 3.0
# Slow DC report cadence is 1 s (Nyquist 0.5 Hz). The FFT swap must cross
# over well below it: centre <= 0.8 x DC Nyquist and the whole raised-cosine
# transition (centre +- DC_TRANSITION_FRAC * centre) below Nyquist. The old
# 1.0 Hz crossover was above the DC Nyquist, i.e. it imported aliased noise.
# DC cadence: MAGIC 1 s (Nyquist 0.5 Hz), LST-1 ~2 s (Nyquist 0.25 Hz).
DC_CUT_HZ = {"m1": 0.3, "m2": 0.3, "lst1": 0.15}
DC_TRANSITION_FRAC = 0.5
# Rolling-median window of the smoothed-DC trend (detrend correction).
DC_DETREND_WINDOW_S = 10.0

@dataclass
class KBOConfig:
    radius_m: float
    distance_au: float
    impact_parameter_m: float = 0.0

@dataclass
class StarConfig:
    temperature_K: float
    angular_radius_mas: float

@dataclass
class BandpassConfig:
    lam_min_nm: float
    lam_max_nm: float
    n_lambda: int

@dataclass
class GridConfig:
    x_max_m: float
    n_x: int

@dataclass
class NumericalConfig:
    n_int: int = 40
    n_r_grid: int = 3000
    n_star_side: int = 32

standard_sampling = 262144 # 2**18
standard_sample_duration = 2 #ns
lightcurve_sampling = 70 # [Hz]. >= 2 * f99 of the most compact star (r_* = 0.01 mas: f99 = 32 Hz at 25 km/s). See README "Signal band"


# ───────────────────────────────────────────────────────────
# Telescope array (single source of truth for the channel <->
# telescope mapping and telescope positions)
# ───────────────────────────────────────────────────────────

@dataclass(frozen=True)
class TelescopeConfig:
    """
    One telescope of the array.

    channel : DAQ channel of the fast (stat-binary) data ("A"/"B"/"C").
    dc_key  : telescope key used by dc_combine.load_dc_series
              ("m1"/"m2"/"lst1").
    east_m, north_m, up_m : local ENU position (m) relative to the
        array reference point (ArrayConfig.latitude_deg / longitude_deg).
    """
    name: str
    channel: str
    dc_key: str
    east_m: float
    north_m: float
    up_m: float = 0.0


@dataclass(frozen=True)
class ArrayConfig:
    """
    The telescope array. Geodetic coordinates anchor the local ENU
    frame the telescope positions are expressed in.
    """
    telescopes: Tuple[TelescopeConfig, ...]
    latitude_deg: float
    longitude_deg: float
    height_m: float

    def by_channel(self, channel: str) -> TelescopeConfig:
        for tel in self.telescopes:
            if tel.channel == channel:
                return tel
        raise KeyError(f"No telescope with channel {channel!r} in this array")

    def channels(self) -> Tuple[str, ...]:
        return tuple(tel.channel for tel in self.telescopes)

    def position_enu_m(self, channel: str) -> np.ndarray:
        tel = self.by_channel(channel)
        return np.array([tel.east_m, tel.north_m, tel.up_m])

    def baseline_enu_m(self, channel_a: str, channel_b: str) -> np.ndarray:
        """ENU vector from telescope channel_a to telescope channel_b."""
        return self.position_enu_m(channel_b) - self.position_enu_m(channel_a)

    def max_baseline_m(self) -> float:
        """Largest pairwise telescope separation (m)."""
        chans = self.channels()
        return max(
            float(np.linalg.norm(self.baseline_enu_m(a, b)))
            for i, a in enumerate(chans) for b in chans[i + 1:]
        )


# MAGIC-1 / MAGIC-2 / LST-1 at the Roque de los Muchachos Observatory.
# Channel mapping (2026 benchmark-star campaign): C = MAGIC-1,
# A = MAGIC-2, B = LST-1.
#
# ENU positions are APPROXIMATE (MAGIC-1 <-> MAGIC-2 ~85 m; LST-1 ~100 m
# from the MAGIC site). At ~25 km/s shadow velocity 100 m is only ~4 ms
# -- small against typical template durations -- so approximate
# positions only pad the coincidence window marginally.
# TODO(user): replace with surveyed telescope coordinates when available.
ORM_ARRAY = ArrayConfig(
    telescopes=(
        TelescopeConfig("MAGIC-1", channel="C", dc_key="m1", east_m=0.0, north_m=0.0),
        TelescopeConfig("MAGIC-2", channel="A", dc_key="m2", east_m=-60.0, north_m=60.0),
        TelescopeConfig("LST-1", channel="B", dc_key="lst1", east_m=80.0, north_m=-70.0),
    ),
    latitude_deg=28.7619,
    longitude_deg=-17.8900,
    height_m=2200.0,
)
