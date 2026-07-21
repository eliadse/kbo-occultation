# kbo_occultation/injection.py
"""
Inject a simulated Fresnel diffraction template into a real light curve.

This is the bridge between the physics simulation (simulation.py) and a
real time series (photometry.LightCurve): it lets you test whether the
detection pipeline (matched_filter.py) can actually recover a known
signal buried in real instrumental noise, rather than only in idealized
noise added to the noiseless simulated curve itself.

Injection is multiplicative (injected = signal * template(t)), matching
the convention already used elsewhere in this package (e.g.
simulation.apply_stellar_disk_2d, where intensity is a dimensionless
transmission fraction that is 1.0 far from the event): whatever real
noise is already present in `signal` is preserved and simply scaled by
the same fractional dip a real occultation would cause.
"""

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .detectability import spatial_to_time


@dataclass
class InjectionResult:
    """
    time, injected_signal : ndarray
        Same shape as the input real light curve.
    template : ndarray
        Per-sample multiplicative factor actually applied (1.0 outside
        the event window), i.e. injected_signal = signal * template.
    t0_s : float
        Time at which the template's x=0 (closest approach) was placed.
    window_s : (float, float)
        Start/end time where `template` deviates from 1.0 by more than a
        small tolerance -- the "truth" window used to score recovery.
    """
    time: np.ndarray
    injected_signal: np.ndarray
    template: np.ndarray
    t0_s: float
    window_s: Tuple[float, float]


def random_injection_time(time: np.ndarray, x_m: np.ndarray, shadow_velocity_mps: float,
                           margin_s: float = 0.0, rng: Optional[np.random.Generator] = None) -> float:
    """
    Pick a random t0_s such that the full template (spanning x_m at
    shadow_velocity_mps) fits inside `time`, with an extra `margin_s` of
    padding on both sides.
    """
    time = np.asarray(time)
    rng = rng or np.random.default_rng()

    half_span_s = float(np.abs(spatial_to_time(x_m, shadow_velocity_mps)).max())
    pad_s = half_span_s + margin_s

    lo, hi = time.min() + pad_s, time.max() - pad_s
    if lo >= hi:
        raise ValueError(
            f"Template span ({2*half_span_s:.3g}s, plus margin) doesn't fit inside "
            f"the light curve duration ({time.max() - time.min():.3g}s)."
        )
    return float(rng.uniform(lo, hi))


def inject_occultation(time: np.ndarray, signal: np.ndarray, x_m: np.ndarray, intensity: np.ndarray,
                        shadow_velocity_mps: float, t0_s: Optional[float] = None,
                        rng: Optional[np.random.Generator] = None, window_tol: float = 1e-3) -> InjectionResult:
    """
    Inject a simulated diffraction pattern (x_m, intensity) into a real
    light curve (time, signal) at time t0_s.

    Parameters
    ----------
    time, signal : array_like
        Real light curve to inject into.
    x_m, intensity : array_like
        Noiseless simulated diffraction pattern, e.g. from
        simulation.compute_lightcurve / OccultationEngine.compute --
        intensity normalised so it is ~1.0 far from the KBO.
    shadow_velocity_mps : float
        KBO shadow velocity across the observer, used to convert x_m to
        time via detectability.spatial_to_time.
    t0_s : float, optional
        Time at which to place the pattern's x=0 (closest approach). If
        None, a random valid time is chosen via random_injection_time
        (pass `rng` for reproducibility).
    window_tol : float
        Fractional deviation from 1.0 above which a sample counts as
        "inside" the event, for the purpose of reporting window_s.

    Returns
    -------
    InjectionResult
    """
    time = np.asarray(time)
    signal = np.asarray(signal)

    if t0_s is None:
        t0_s = random_injection_time(time, x_m, shadow_velocity_mps, rng=rng)

    t_template = spatial_to_time(x_m, shadow_velocity_mps) + t0_s
    # Interpolate in float64, then match the data's dtype so a float32
    # light curve stays float32 downstream (halves per-trial memory).
    template = np.interp(time, t_template, intensity, left=1.0, right=1.0)
    template = template.astype(signal.dtype, copy=False)

    injected_signal = signal * template

    inside = np.abs(template - 1.0) > window_tol
    if np.any(inside):
        idx = np.where(inside)[0]
        window_s = (float(time[idx[0]]), float(time[idx[-1]]))
    else:
        window_s = (t0_s, t0_s)

    return InjectionResult(
        time=time,
        injected_signal=injected_signal,
        template=template,
        t0_s=t0_s,
        window_s=window_s,
    )
