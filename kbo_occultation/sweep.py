# kbo_occultation/sweep.py
"""
Parameter-sweep utilities for testing how KBO size, KBO distance, impact
parameter, stellar angular size, and instrument bandpass/filter affect the
occultation diffraction signal.

Design note
-----------
For a fixed (KBO radius, KBO distance, bandpass/response, star
temperature), the per-wavelength Fresnel radial intensity profile does
NOT depend on impact parameter or stellar angular size -- those two only
affect the stellar-disk convolution step that happens afterwards.
OccultationEngine (in simulation.py) caches the radial profile keyed on
(radius, distance) on a fixed radial grid, so the efficient way to sweep
is:

    1. Group by (bandpass/response, star temperature) -- these determine
       the wavelength grid and spectral weights. Build ONE engine per
       group.
    2. Within a group, loop over KBO radius/distance (each combination
       computes and caches a new Fresnel profile), and impact parameter /
       stellar size (each of these reuses the cached profile and only
       redoes the cheap stellar-disk convolution).

`run_parameter_sweep` below does this grouping for you.
"""

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from .config import KBOConfig, StarConfig, BandpassConfig, GridConfig, NumericalConfig
from .physics import AU_m, mas_to_rad
from .simulation import OccultationEngine


@dataclass
class SweepPoint:
    """
    One evaluated light curve, tagged with the parameters that produced
    it, plus a couple of cheap summary metrics computed up front (so
    they're available even if you drop the full curves to save memory).
    """
    radius_m: float
    distance_au: float
    impact_parameter_m: float
    star_temperature_K: float
    star_angular_radius_mas: float
    band_label: str
    depth: float
    min_intensity: float
    x_m: Optional[np.ndarray] = None
    intensity: Optional[np.ndarray] = None

    def summary(self) -> dict:
        """Dict of everything except the (possibly large) curve arrays."""
        return {
            "radius_m": self.radius_m,
            "distance_au": self.distance_au,
            "impact_parameter_m": self.impact_parameter_m,
            "star_temperature_K": self.star_temperature_K,
            "star_angular_radius_mas": self.star_angular_radius_mas,
            "band": self.band_label,
            "depth": self.depth,
            "min_intensity": self.min_intensity,
        }


def _band_label(band: BandpassConfig) -> str:
    return f"{band.lam_min_nm:g}-{band.lam_max_nm:g}nm"


def _required_r_max_m(x_max_m, star_angular_radii_mas, distances_au, impact_parameters_m):
    """
    Upper bound on the radial grid needed to cover every combination of
    star size, distance, and impact parameter in a sweep, following the
    same geometry OccultationEngine.compute() uses internally.
    """
    r_star_max_m = max(star_angular_radii_mas) * mas_to_rad * max(distances_au) * AU_m
    b_max_m = max(abs(b) for b in impact_parameters_m)
    return np.sqrt((x_max_m + r_star_max_m)**2 + (r_star_max_m + b_max_m)**2)


def run_parameter_sweep(
    radii_m: Sequence[float],
    distances_au: Sequence[float],
    impact_parameters_m: Sequence[float],
    star_angular_radii_mas: Sequence[float],
    star_temperature_K: float,
    bands: Sequence[BandpassConfig],
    grid: GridConfig,
    numerics: NumericalConfig,
    responses: Optional[Dict[str, Callable]] = None,
    r_max_m: Optional[float] = None,
    weight_rel_tol: float = 1e-6,
    keep_curves: bool = False,
) -> List[SweepPoint]:
    """
    Sweep the diffraction signal over KBO radius, KBO distance, impact
    parameter, stellar angular size, and bandpass/filter.

    Parameters
    ----------
    radii_m, distances_au, impact_parameters_m, star_angular_radii_mas :
        Sequences of values to sweep over.
    star_temperature_K : float
        Held fixed across the sweep (spectral weights depend on it; add
        an outer loop yourself if you also want to sweep temperature).
    bands : sequence of BandpassConfig
        One or more bandpasses/filters to test.
    responses : dict, optional
        Maps a band's label (``f"{lam_min_nm:g}-{lam_max_nm:g}nm"``) to a
        response callable R(lambda_nm), e.g. built with
        instruments.build_response_function / combine_responses. If a
        band has no entry here, its own trapezoidal filter_transmission
        is used instead.
    r_max_m : float, optional
        Fixed outer radius of the radial evaluation grid, shared by every
        point in the sweep -- must be large enough to cover x_max plus
        the largest impact parameter and stellar radius involved. If
        omitted, it's computed automatically from the other arguments.
        Passing it explicitly (a bit larger than the auto value) avoids
        any mid-sweep grid-extension warnings if you later add points.
    weight_rel_tol : float
        See OccultationEngine.compute.
    keep_curves : bool
        If False, x_m/intensity are dropped after each point's summary
        metrics (depth, min_intensity) are computed, to save memory on
        large sweeps.

    Returns
    -------
    list of SweepPoint
    """
    responses = responses or {}

    if r_max_m is None:
        r_max_m = _required_r_max_m(grid.x_max_m, star_angular_radii_mas, distances_au, impact_parameters_m)

    results: List[SweepPoint] = []

    for band in bands:
        label = _band_label(band)
        response = responses.get(label)

        star_for_engine = StarConfig(star_temperature_K,star_angular_radii_mas[0])
        engine = OccultationEngine(star_for_engine, band, grid, numerics, response, r_max_m)

        for distance_au in distances_au:
            for radius_m in radii_m:
                for impact_m in impact_parameters_m:
                    kbo = KBOConfig(radius_m, distance_au, impact_m)
                    for star_size_mas in star_angular_radii_mas:
                        x_m, I = engine.compute(kbo, star_size_mas, weight_rel_tol)
                        depth = 1.0 - float(np.min(I))
                        point = SweepPoint(
                            radius_m=radius_m,
                            distance_au=distance_au,
                            impact_parameter_m=impact_m,
                            star_temperature_K=star_temperature_K,
                            star_angular_radius_mas=star_size_mas,
                            band_label=label,
                            depth=depth,
                            min_intensity=float(np.min(I)),
                            x_m=x_m.copy() if keep_curves else None,
                            intensity=I.copy() if keep_curves else None,
                        )
                        results.append(point)

    return results


def to_records(points: Sequence[SweepPoint]) -> List[dict]:
    """Plain list-of-dicts summary, no pandas required."""
    return [p.summary() for p in points]


def to_dataframe(points: Sequence[SweepPoint]):
    """
    Convert sweep results to a pandas DataFrame of summary metrics
    (one row per parameter combination; curve arrays are not included).
    Requires pandas to be installed.
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "to_dataframe() requires pandas (pip install pandas); "
            "use to_records() for a plain list of dicts instead."
        ) from e
    return pd.DataFrame(to_records(points))
