# kbo_occultation/analysis/spectra.py
"""
Noise spectra of light curves, gap-aware, and the two noise numbers that
matter for a matched-filter search:

- inband_rms: noise rms inside a frequency band (the diffraction-signal
  band from analysis.bands);
- matched_filter_sigma: sqrt(int |T|^2 S df / int |T|^2 df), the noise a
  filter with template T actually sees in a noise PSD S. It equals the
  empirical scatter of the matched-filter output
  (matched_filter.empirical_filter_sigma) for stationary Gaussian noise,
  and is what a white-noise sigma / sqrt(sum T^2) hides when the noise is
  red.

All spectra are of the fractional flux (signal / mean - 1), one-sided, in
1/Hz.
"""

from typing import Dict, Optional, Tuple

import numpy as np
from scipy.signal import welch

SEGMENT_SAMPLES = 2 ** 16          # ~34 s at 0.524 ms: 0.03 Hz resolution


def fractional_flux(signal: np.ndarray) -> np.ndarray:
    signal = np.asarray(signal, dtype=float)
    return signal / np.nanmean(signal) - 1.0


def segmented_welch(time, signal, fs, nperseg_target: int = SEGMENT_SAMPLES,
                    min_run_samples: int = 1024, return_info: bool = False):
    """
    Welch PSD that does not straddle gaps: the series is split wherever
    the time step exceeds 1.5 / fs (e.g. after despike with fill="remove"),
    each run gets its own Welch estimate, and the estimates are averaged
    weighted by run length, on the common grid of ``nperseg_target``.

    Runs shorter than ``min_run_samples`` are skipped. Runs shorter than
    ``nperseg_target`` use a shorter segment (coarser frequency
    resolution, interpolated onto the common grid); this is reported in
    ``info`` so the caller can flag it.

    Returns (freq, psd), or (freq, psd, info) with ``return_info=True``;
    info = {"nperseg": longest segment actually used, "short": True if any
    run used less than nperseg_target, "n_runs", "used_fraction"}.
    """
    time = np.asarray(time)
    signal = np.asarray(signal)
    dt = 1.0 / fs
    break_idx = np.where(np.diff(time) > 1.5 * dt)[0] + 1
    bounds = np.concatenate([[0], break_idx, [len(time)]])

    freq_common = np.linspace(0, fs / 2, nperseg_target // 2 + 1)
    accum, weight_total, used, longest, n_runs = None, 0, 0, 0, 0
    for start, stop in zip(bounds[:-1], bounds[1:]):
        run = signal[start:stop]
        if len(run) < min_run_samples:
            continue
        nperseg = min(nperseg_target, len(run))
        f, pxx = welch(run, fs=fs, window="hann", nperseg=nperseg,
                       noverlap=nperseg // 2, scaling="density")
        pxx = np.interp(freq_common, f, pxx)
        accum = pxx * len(run) if accum is None else accum + pxx * len(run)
        weight_total += len(run)
        used += len(run)
        longest = max(longest, nperseg)
        n_runs += 1
    if accum is None:
        raise ValueError("no run long enough for a Welch estimate")
    psd = accum / weight_total
    if not return_info:
        return freq_common, psd
    info = {"nperseg": longest, "short": longest < nperseg_target or n_runs > 1 and any(
                stop - start < nperseg_target for start, stop in zip(bounds[:-1], bounds[1:])
                if stop - start >= min_run_samples),
            "n_runs": n_runs, "used_fraction": used / len(time)}
    return freq_common, psd, info


def blockwise_welch(time, signal, fs, n_block: int, nperseg: int = 1024, average: str = "median",
                    min_coverage: float = 0.5):
    """
    Welch PSD of ``n_block`` equal time blocks of one series, gap-aware as ``segmented_welch``
    (a block's runs are split where the time step exceeds 1.5 / fs, runs shorter than ``nperseg``
    are skipped, the run estimates are averaged weighted by length). ``average="median"`` (the
    default) makes each run's estimate robust to noise bursts, like the median-based scatter of the
    matched-filter normalisation.

    Returns (freq, psd): ``psd`` is (n_block, nperseg // 2 + 1); a block that has less than
    ``min_coverage`` of its samples in usable runs is NaN.
    """
    time, signal = np.asarray(time), np.asarray(signal)
    dt = 1.0 / fs
    edges = np.linspace(time[0], time[-1] + dt, n_block + 1)
    bounds = np.searchsorted(time, edges)
    freq = np.linspace(0.0, fs / 2.0, nperseg // 2 + 1)
    out = np.full((n_block, len(freq)), np.nan)
    for k in range(n_block):
        a, b = bounds[k], bounds[k + 1]
        if b - a < nperseg:
            continue
        cuts = np.concatenate([[0], np.where(np.diff(time[a:b]) > 1.5 * dt)[0] + 1, [b - a]])
        accum, used = 0.0, 0
        for s, e in zip(cuts[:-1], cuts[1:]):
            if e - s < nperseg:
                continue
            _, pxx = welch(signal[a + s:a + e], fs=fs, window="hann", nperseg=nperseg,
                           noverlap=nperseg // 2, scaling="density", average=average)
            accum = accum + pxx * (e - s)
            used += e - s
        expected = (edges[k + 1] - edges[k]) * fs
        if used and used >= min_coverage * expected:
            out[k] = accum / used
    return freq, out


def lightcurve_psd(lc, nperseg_target: int = SEGMENT_SAMPLES, return_info: bool = False):
    """segmented_welch of a LightCurve's fractional flux."""
    return segmented_welch(lc.time, fractional_flux(lc.signal), 1.0 / lc.dt,
                           nperseg_target, return_info=return_info)


def _trapz(y, x):
    return float((getattr(np, "trapezoid", None) or np.trapz)(y, x))


def inband_rms(lc, band: Tuple[float, float], nperseg_target: int = SEGMENT_SAMPLES,
               return_info: bool = False):
    """
    RMS of the fractional flux inside ``band`` = (f_lo, f_hi) Hz, from the
    gap-aware Welch PSD. With ``return_info`` also returns the Welch info
    dict (``short`` flags runs too short for ``nperseg_target``).
    """
    f, psd, info = lightcurve_psd(lc, nperseg_target, return_info=True)
    sel = (f >= band[0]) & (f <= band[1])
    rms = float(np.sqrt(_trapz(psd[sel], f[sel]))) if sel.sum() > 1 else float("nan")
    return (rms, info) if return_info else rms


def matched_filter_sigma(freq: np.ndarray, psd: np.ndarray, template_delta: np.ndarray,
                         fs: float) -> float:
    """
    Noise seen by a matched filter: sqrt( int |T|^2 S df / int |T|^2 df ).

    Parameters
    ----------
    freq, psd : the PSD (one-sided, fractional-flux^2 / Hz) on the grid of
        ``segmented_welch`` (spacing fs / nperseg).
    template_delta : signal variation of the template, template - 1, sampled at 1 / fs.

    Same units as the per-sample sigma that gives the same SNR for the
    unit-norm filter (matched_filter.MatchedFilterResult.amp_sigma). For
    white noise of per-sample sigma s, S = 2 s^2 / fs and the result is s.
    """
    n_fft = 2 * (len(freq) - 1)
    T = np.fft.rfft(np.asarray(template_delta, dtype=float), n=n_fft)
    w = np.abs(T) ** 2
    return float(np.sqrt(np.sum(w * psd) * (fs / 2.0) / np.sum(w)))


def noise_band_numbers(lc, band: Tuple[float, float], template_delta: Optional[np.ndarray] = None,
                       nperseg_target: int = SEGMENT_SAMPLES) -> Dict[str, float]:
    """
    One Welch estimate -> {"inband_rms", "welch_short"} and, when a template
    variation (template - 1, sampled at lc.dt) is given, "mf_sigma". All in
    fractional-flux units; "welch_short" is 1.0 when some run was too short
    for ``nperseg_target`` (coarser frequency resolution than intended).
    """
    f, psd, info = lightcurve_psd(lc, nperseg_target, return_info=True)
    sel = (f >= band[0]) & (f <= band[1])
    out = {"inband_rms": float(np.sqrt(_trapz(psd[sel], f[sel]))) if sel.sum() > 1 else float("nan"),
           "welch_short": float(info["short"])}
    if template_delta is not None:
        out["mf_sigma"] = matched_filter_sigma(f, psd, template_delta, 1.0 / lc.dt)
    return out
