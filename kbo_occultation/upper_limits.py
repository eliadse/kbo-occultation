# kbo_occultation/upper_limits.py
"""
Upper limits on the cumulative surface density of KBOs from a null
occultation search.

Formalism (TAOS-style; Nihei et al. 2007, Zhang et al. 2008/2013): each
star observed for a live time T sweeps out a solid angle of "shadow
opportunity" per unit radius class,

    Omega(r) = [H(r) / D] * [v_rel * T / D]      [steradians],

where H(r) is the effective shadow width (diffraction-aware event
cross-section), D the occulter distance, and v_rel the shadow velocity.
The expected number of detected events for a cumulative surface density
N(>r0) [deg^-2] of objects with radius > r0 is at least

    mu >= N(>r0) * sum_i eps_i(r0) * Omega_i(r0) * (180/pi)^2,

because both the detection efficiency eps(r) and H(r) are
non-decreasing in r -- evaluating everything at r0 (as if every object
had exactly radius r0) is the standard conservative treatment for a
cumulative limit. With n observed events, the Poisson upper limit
mu_UL(n, CL) then bounds

    N(>r0) < mu_UL / sum_i [ eps_i(r0) * Omega_i(r0) * (180/pi)^2 ].

For zero detections at 95% CL, mu_UL = -ln(0.05) = 2.996.
"""

from dataclasses import dataclass, field
from typing import Dict, Sequence

import numpy as np
from scipy.stats import chi2

AU_M = 1.495978707e11


@dataclass
class EfficiencyCurve:
    """
    Detection efficiency vs radius, e.g. from
    injection_batch.efficiency_curve. Calling it interpolates linearly,
    clips to [0, 1], and returns 0 outside the measured radius range
    (no extrapolated sensitivity is claimed).
    """
    radii_m: np.ndarray
    efficiency: np.ndarray

    def __call__(self, radius_m):
        radius_m = np.asarray(radius_m, dtype=float)
        eps = np.interp(radius_m, self.radii_m, self.efficiency, left=0.0, right=0.0)
        return np.clip(eps, 0.0, 1.0)


@dataclass
class SurveyExposure:
    """
    One star's contribution to the survey exposure.

    live_time_s : effective observing time after dead-time/despike
        removal.
    shadow_velocity_mps : the velocity the efficiency curve was
        measured at / the survey is bookkept at (typically the
        opposition value, kinematics.shadow_velocity_opposition).
    star_radius_projected_m : stellar angular radius projected at the
        occulter distance (angular_radius_rad * D); broadens the
        effective shadow width.
    wavelength_m : band-representative wavelength for the Fresnel scale.
    """
    star_name: str
    live_time_s: float
    shadow_velocity_mps: float
    distance_au: float
    efficiency: EfficiencyCurve
    star_radius_projected_m: float = 0.0
    wavelength_m: float = 500e-9


@dataclass
class UpperLimitResult:
    radii_m: np.ndarray
    n_upper_deg2: np.ndarray
    mu_upper: float
    confidence: float
    per_star_omega_deg2: Dict[str, np.ndarray] = field(default_factory=dict)


def fresnel_scale_m(wavelength_m: float, distance_au: float) -> float:
    """F = sqrt(lambda * D / 2). ~1.3 km at 500 nm, 43 AU."""
    return float(np.sqrt(wavelength_m * distance_au * AU_M / 2.0))


def effective_shadow_width_m(radius_m, distance_au: float, wavelength_m: float,
                             star_radius_projected_m: float = 0.0):
    """
    Diffraction-aware event cross-section (Nihei et al. 2007; the TAOS
    convention):

        H(r) = 2 * ((sqrt(3) * F)**1.5 + r**1.5)**(2/3) + 2 * r_star

    Limits: H -> 2*sqrt(3)*F + 2*r_star for r << F (diffraction-
    dominated), H -> 2*r for r >> F (geometric).
    """
    radius_m = np.asarray(radius_m, dtype=float)
    F = fresnel_scale_m(wavelength_m, distance_au)
    return (2.0 * ((np.sqrt(3.0) * F) ** 1.5 + radius_m ** 1.5) ** (2.0 / 3.0)
            + 2.0 * star_radius_projected_m)


def poisson_upper_limit(n_observed: int = 0, confidence: float = 0.95) -> float:
    """
    Classical Poisson upper limit on the mean, mu_UL such that
    P(<= n_observed | mu_UL) = 1 - confidence:

        mu_UL = 0.5 * chi2.ppf(confidence, 2 * (n_observed + 1))

    = 2.996 for zero events at 95% CL.
    """
    return float(0.5 * chi2.ppf(confidence, 2 * (n_observed + 1)))


def exposure_omega_deg2(exposure: SurveyExposure, radii_m: np.ndarray) -> np.ndarray:
    """
    Solid angle swept by detectable shadows for one star, per radius,
    in deg^2 (before the efficiency factor).
    """
    radii_m = np.asarray(radii_m, dtype=float)
    D_m = exposure.distance_au * AU_M
    H = effective_shadow_width_m(radii_m, exposure.distance_au, exposure.wavelength_m,
                                 exposure.star_radius_projected_m)
    omega_sr = (H / D_m) * (exposure.shadow_velocity_mps * exposure.live_time_s / D_m)
    return omega_sr * (180.0 / np.pi) ** 2


def cumulative_density_upper_limit(exposures: Sequence[SurveyExposure],
                                   radii_m: np.ndarray,
                                   confidence: float = 0.95,
                                   n_observed: int = 0) -> UpperLimitResult:
    """
    Upper limit on the cumulative surface density N(> r0) [deg^-2] as a
    function of the threshold radius r0, from a set of per-star
    exposures with measured efficiency curves.
    """
    radii_m = np.asarray(radii_m, dtype=float)
    mu_upper = poisson_upper_limit(n_observed, confidence)

    per_star_omega: Dict[str, np.ndarray] = {}
    effective_omega = np.zeros_like(radii_m)
    for exp in exposures:
        omega = exposure_omega_deg2(exp, radii_m)
        per_star_omega[exp.star_name] = omega
        effective_omega += exp.efficiency(radii_m) * omega

    with np.errstate(divide="ignore"):
        n_upper = np.where(effective_omega > 0, mu_upper / effective_omega, np.inf)

    return UpperLimitResult(
        radii_m=radii_m,
        n_upper_deg2=n_upper,
        mu_upper=mu_upper,
        confidence=confidence,
        per_star_omega_deg2=per_star_omega,
    )
