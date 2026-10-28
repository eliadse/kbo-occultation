# kbo_occultation/analysis/surrogates.py
"""
Gaussian surrogates of a real light curve: same noise spectrum, same
gaps, none of the non-Gaussian structure (gusts, glitches, transients).

A search on surrogates gives the false-alarm rate that stationary
Gaussian noise with the *measured* spectrum produces. A real-data rate
above it points to transient, non-Gaussian mimics.
"""

from typing import Iterator, Optional

import numpy as np

from ..observations.photometry import LightCurve


def _runs(time: np.ndarray, dt: float):
    breaks = np.where(np.diff(time) > 1.5 * dt)[0] + 1
    bounds = np.concatenate([[0], breaks, [len(time)]])
    return list(zip(bounds[:-1], bounds[1:]))


def phase_randomize_array(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """
    Keep the amplitude spectrum of ``x`` and draw uniform random phases
    (the DC bin and, for even length, Nyquist bin keep a real value). The
    mean and the periodogram are preserved exactly; the output is a
    stationary Gaussian-like series.
    """
    x = np.asarray(x, dtype=float)
    spec = np.fft.rfft(x - x.mean())
    phases = np.exp(2j * np.pi * rng.random(len(spec)))
    phases[0] = 1.0
    if len(x) % 2 == 0:
        phases[-1] = 1.0
    return np.fft.irfft(np.abs(spec) * phases, n=len(x)) + x.mean()


def phase_randomized_lightcurve(lc: LightCurve, rng: Optional[np.random.Generator] = None) -> LightCurve:
    """
    One surrogate of ``lc``. Each contiguous run (split at time gaps, e.g.
    the holes left by despike with fill="remove") is phase-randomised on
    its own, so the gap mask and every run's length and mean are kept
    exactly and no spectral content leaks across a gap. The time axis is
    shared with ``lc``.
    """
    rng = rng or np.random.default_rng()
    signal = np.asarray(lc.signal, dtype=float)
    out = signal.copy()
    for start, stop in _runs(np.asarray(lc.time), lc.dt):
        if stop - start >= 8:
            out[start:stop] = phase_randomize_array(signal[start:stop], rng)
    return LightCurve(lc.time, out, meta=dict(lc.meta, surrogate=True))


def surrogate_lightcurves(lc: LightCurve, n: int, seed: Optional[int] = None) -> Iterator[LightCurve]:
    """``n`` independent surrogates of ``lc`` (reproducible with ``seed``)."""
    rng = np.random.default_rng(seed)
    for _ in range(n):
        yield phase_randomized_lightcurve(lc, rng)
