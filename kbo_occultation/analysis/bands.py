# kbo_occultation/analysis/bands.py
"""
Frequency band occupied by the diffraction signal.

The band is described by several numbers rather than a single
"Nyquist" threshold, because any single cumulative-power cut-off is
arbitrary:

f05, f50, f90, f95, f99 : temporal frequencies (Hz) below which that
    fraction of the power of the signal (1 - I) lies.
f_low : f05. A high-pass with cutoff below f_low removes < 5 % of the
    signal power.
f_mf_hi : f95, the upper edge of the band holding the central 90 %
    ([f05, f95]) of the matched-filter weight int |T(f)|^2 df. For white
    noise the matched-filter SNR^2 is proportional to this integral, so
    it is the band the filter is sensitive in. The width is f_mf_hi - f_low.
chord_time_s : H / v, the event duration from effective_shadow_width_m
    (Nihei et al. 2007 diffraction-aware width plus the stellar disk).

All spectra are of (1 - I) on the native spatial grid, converted to
Hz with f = f_x * v.
"""

from functools import lru_cache
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from ..config import (BandpassConfig, GridConfig, KBOConfig, NumericalConfig,
                      StarConfig)
from ..simulation.detectability import (cumulative_quantile_frequency,
                                        spatial_power_spectrum)
from ..simulation.engine import OccultationEngine
from ..simulation.physics import AU_m, mas_to_rad
from .upper_limits import effective_shadow_width_m

DEFAULT_QUANTILES = (0.05, 0.50, 0.90, 0.95, 0.99)


def signal_band(star_mas: float, radii_m: Sequence[float], band: BandpassConfig,
                distance_au: float, v_mps: float,
                quantiles: Sequence[float] = DEFAULT_QUANTILES,
                impact_parameter_m: float = 0.0, temperature_K: float = 10000.0,
                response: Optional[Callable] = None,
                grid: Optional[GridConfig] = None,
                numerics: Optional[NumericalConfig] = None,
                interpolate: bool = True, pad_factor: int = 4) -> List[Dict[str, float]]:
    """
    Band descriptors of the diffraction signal for one stellar size, one
    distance and one velocity, for every radius in ``radii_m``.

    Returns one dict per radius with keys ``radius_m``, ``f<NN>_hz`` for
    each quantile (e.g. ``f99_hz``), ``f_low_hz``, ``f_mf_hi_hz``,
    ``bandwidth_hz``, ``chord_time_s`` and ``depth``.
    """
    grid = grid or GridConfig(x_max_m=30000.0, n_x=12000)
    numerics = numerics or NumericalConfig(n_r_grid=6000, n_star_side=128)

    r_star_m = star_mas * mas_to_rad * distance_au * AU_m
    r_max = 1.2 * np.sqrt((grid.x_max_m + r_star_m) ** 2
                          + (r_star_m + abs(impact_parameter_m)) ** 2)
    engine = OccultationEngine(StarConfig(temperature_K, star_mas), band, grid, numerics,
                               response=response, r_max_m=r_max)
    lam_eff_m = float(np.sum(engine.lambdas_nm * engine.weights)) * 1e-9

    rows = []
    for radius in radii_m:
        kbo = KBOConfig(float(radius), distance_au, impact_parameter_m)
        x_m, intensity = engine.compute(kbo, star_mas)
        fx, power = spatial_power_spectrum(x_m, intensity, baseline=1.0, pad_factor=pad_factor)
        row = {"radius_m": float(radius), "star_mas": star_mas, "distance_au": distance_au,
               "v_mps": v_mps, "depth": float(1.0 - intensity.min())}
        for q in quantiles:
            row[f"f{int(round(100 * q)):02d}_hz"] = (
                cumulative_quantile_frequency(fx, power, q, interpolate) * v_mps)
        row["f_low_hz"] = row.get("f05_hz", np.nan)
        row["f_mf_hi_hz"] = row.get("f95_hz", np.nan)
        row["bandwidth_hz"] = row["f_mf_hi_hz"] - row["f_low_hz"]
        row["chord_time_s"] = float(effective_shadow_width_m(radius, distance_au, lam_eff_m,
                                                             r_star_m) / v_mps)
        rows.append(row)
    return rows


# Band (Hz) used for the in-band noise numbers of the SNR report. Lower edge:
# ~f05 of the compact-star (r* = 0.03 mas) pattern; upper edge: f95-f99 range
# of r* = 0.01-0.03 mas. See README "Signal band".
SIGNAL_BAND_HZ = (0.2, 25.0)


@lru_cache(maxsize=8)
def default_template_residual(dt_s: float, radius_m: float = 100.0, star_mas: float = 0.03,
                              distance_au: float = 40.0, v_mps: float = 25e3) -> np.ndarray:
    """
    Reference template (I - 1) sampled at ``dt_s``: a 100 m KBO at 40 AU
    crossing a 0.03 mas star at 25 km/s, b = 0, in the 320-450 nm band.
    Light-weight numerics (the shape only sets the weighting of the noise
    PSD). Cached per ``dt_s``.
    """
    half = 4000.0
    grid = GridConfig(half, int(2 * half / 5.0))
    engine = OccultationEngine(StarConfig(10000.0, star_mas), BandpassConfig(320, 450, 21), grid,
                               NumericalConfig(n_r_grid=2000, n_star_side=32),
                               r_max_m=1.3 * (half + star_mas * mas_to_rad * distance_au * AU_m) * 1.5)
    x_m, intensity = engine.compute(KBOConfig(radius_m, distance_au, 0.0), star_mas)
    t = x_m / v_mps
    grid_t = np.arange(t.min(), t.max(), dt_s)
    return np.interp(grid_t, t, intensity) - 1.0
