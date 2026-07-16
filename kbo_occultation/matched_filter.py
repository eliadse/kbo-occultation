# kbo_occultation/matched_filter.py
"""
Blind matched-filter search for an occultation-shaped dip in a light
curve, at an unknown location.

detectability.matched_filter_snr computes the SNR of a dip *already
aligned* with a template -- useful for mapping out detectability across
a parameter grid where you control both. This module does the harder,
more realistic thing: slide the template across the whole series and
report an SNR(t) curve, the way a real search over survey data would,
with no assumption about where (or whether) an event is.
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import correlate, find_peaks


@dataclass
class MatchedFilterResult:
    time: np.ndarray
    snr: np.ndarray
    sigma: float
    baseline: np.ndarray
    dt_s: float
    template_duration_s: float
    # Kept around so shape_veto_chi2 can re-fit the template against any
    # candidate's native-resolution data window without redoing the
    # detrending/resampling step.
    data_time: np.ndarray
    residual_data: np.ndarray
    residual_template: np.ndarray


@dataclass
class Candidate:
    time_s: float
    snr: float
    chi2_reduced: float


def robust_sigma(x: np.ndarray) -> float:
    """
    MAD-based noise sigma, robust to a brief injected/real dip (or a
    handful of outliers) not biasing the estimate the way a plain
    ``std`` would over a long series with a short event.
    """
    x = np.asarray(x)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def compute_baseline_and_sigma(time: np.ndarray, signal: np.ndarray, template_duration_s: float,
                                detrend_window_s: Optional[float] = None) -> Tuple[np.ndarray, float]:
    """
    Rolling-median baseline (scipy.ndimage.median_filter) and the
    resulting residual's robust_sigma, factored out of
    sliding_matched_filter_snr so it can be computed once and reused
    (e.g. injection_batch.py caches this per noise-reduction method
    instead of recomputing it -- by far the most expensive step here --
    for every trial).

    Parameters
    ----------
    time, signal : array_like
        Light curve, on a uniform time grid.
    template_duration_s : float
        Used only to set the default detrend_window_s (20x this).
    detrend_window_s : float, optional
        Width of the rolling-median baseline. Defaults to 20x
        template_duration_s, wide enough to not distort the event
        itself.

    Returns
    -------
    baseline : ndarray, same shape as signal.
    sigma : float, robust_sigma of (signal/baseline - 1).
    """
    time = np.asarray(time)
    signal = np.asarray(signal)
    dt = float(np.median(np.diff(time)))

    if detrend_window_s is None:
        detrend_window_s = 20.0 * template_duration_s
    M = int(round(template_duration_s / dt))
    window_samples = max(int(round(detrend_window_s / dt)), M)
    baseline = median_filter(signal, size=window_samples, mode="reflect")
    sigma = robust_sigma(signal / baseline - 1.0)

    return baseline, sigma


def sliding_matched_filter_snr(time: np.ndarray, signal: np.ndarray,
                                template_time: np.ndarray, template_intensity: np.ndarray,
                                detrend_window_s: Optional[float] = None,
                                baseline: Optional[np.ndarray] = None,
                                sigma: Optional[float] = None) -> MatchedFilterResult:
    """
    Slide (template_time, template_intensity) -- a noiseless simulated
    diffraction pattern already converted to time, e.g. via
    detectability.spatial_to_time -- across (time, signal) and compute
    the matched-filter SNR at every lag.

    Parameters
    ----------
    time, signal : array_like
        Light curve to search (real or injected-real), on a uniform
        time grid.
    template_time, template_intensity : array_like
        Noiseless template, time-relative (need not be centered on
        t=0; only its shape and native resolution matter here).
    detrend_window_s : float, optional
        Width of the rolling-median baseline removed before matching,
        to strip slow drifts (elevation, clouds) without needing to
        know the event scale ahead of time. Defaults to 20x the
        template's duration, wide enough to not distort the event
        itself. Ignored if `baseline`/`sigma` are both given.
    baseline, sigma : optional
        Precomputed compute_baseline_and_sigma(...) output, to skip
        recomputing the rolling-median baseline (the expensive step --
        the FFT correlate and per-candidate chi2 veto are cheap by
        comparison). `baseline` must be the same length as `signal`.
        If omitted (the default), both are computed fresh here exactly
        as before -- every existing caller is unaffected.

    Returns
    -------
    MatchedFilterResult
    """
    time = np.asarray(time)
    signal = np.asarray(signal)
    template_time = np.asarray(template_time)
    template_intensity = np.asarray(template_intensity)

    dt = float(np.median(np.diff(time)))

    # --- resample the template onto the data's native cadence ---
    template_duration_s = float(template_time.max() - template_time.min())
    template_grid = np.arange(template_time.min(), template_time.max(), dt)
    template_resampled = np.interp(template_grid, template_time, template_intensity)
    M = len(template_resampled)

    # --- detrend: rolling-median baseline, then normalise ---
    if baseline is None or sigma is None:
        baseline, sigma = compute_baseline_and_sigma(time, signal, template_duration_s, detrend_window_s)
    norm_signal = signal / baseline

    # --- matched filter: FFT cross-correlation of the residuals ---
    residual_data = norm_signal - 1.0
    residual_template = template_resampled - 1.0
    template_norm = float(np.sqrt(np.sum(residual_template**2)))

    correlation = correlate(residual_data, residual_template, mode="valid", method="fft")

    snr = correlation / (template_norm * sigma)

    # correlation[i] is the match for the window signal[i:i+M]; tag it
    # at that window's center time.
    center_offset = M // 2
    snr_time = time[center_offset:center_offset + len(snr)]

    return MatchedFilterResult(
        time=snr_time,
        snr=snr,
        sigma=sigma,
        baseline=baseline,
        dt_s=dt,
        template_duration_s=template_duration_s,
        data_time=time,
        residual_data=residual_data,
        residual_template=residual_template,
    )


def shape_veto_chi2(result: MatchedFilterResult, candidate_time_s: float) -> float:
    """
    Reduced chi-squared of the data around candidate_time_s against the
    *whole* template shape (fringes included), after fitting only its
    amplitude by least squares.

    The plain SNR from sliding_matched_filter_snr is just a correlation:
    it responds to anything that overlaps the template's gross
    depth/duration, including an unrelated dip that happens to share
    that envelope but has none of the diffraction pattern's fine
    structure (secondary fringes either side of the main minimum). This
    chi-squared veto (the same idea as the CBC chi^2 veto used in
    gravitational-wave matched-filter searches, Allen 2005) checks the
    *residuals* after best-fitting the amplitude: a genuine event has to
    reproduce the fringes too, so its residuals stay noise-like and
    chi2/dof stays near 1; a dip that only mimics the envelope leaves
    the unmatched fringe structure behind as excess residual power, and
    chi2/dof rises well above 1.

    Parameters
    ----------
    result : MatchedFilterResult
        Output of sliding_matched_filter_snr on the same data.
    candidate_time_s : float
        Time to test, e.g. from a Candidate returned by find_candidates.

    Returns
    -------
    float
        chi2 / dof, dof = (number of samples in the template window - 1).
    """
    M = len(result.residual_template)
    center_idx = int(np.searchsorted(result.data_time, candidate_time_s))
    start = max(0, min(center_idx - M // 2, len(result.residual_data) - M))
    window = result.residual_data[start:start + M]

    template_norm_sq = float(np.sum(result.residual_template**2))
    amplitude = float(np.sum(window * result.residual_template) / template_norm_sq)

    residuals = window - amplitude * result.residual_template
    dof = M - 1
    return float(np.sum(residuals**2) / result.sigma**2 / dof)


def find_candidates(result: MatchedFilterResult, snr_threshold: float = 5.0,
                     min_separation_s: Optional[float] = None,
                     max_chi2_reduced: Optional[float] = None) -> List[Candidate]:
    """
    Peaks of result.snr above snr_threshold, separated by at least
    min_separation_s (defaults to the template duration, so sidelobes
    of the same event don't count as separate candidates), sorted by
    descending SNR.

    Each candidate is also scored with shape_veto_chi2 (see its
    docstring). Pass max_chi2_reduced to drop candidates whose shape
    doesn't actually match the template's fine structure -- e.g. real
    but unrelated dips that only share the template's gross
    envelope/duration -- even though their correlation SNR cleared
    snr_threshold.
    """
    if min_separation_s is None:
        min_separation_s = result.template_duration_s

    distance = max(int(round(min_separation_s / result.dt_s)), 1)
    peak_idx, _ = find_peaks(result.snr, height=snr_threshold, distance=distance)

    candidates = [
        Candidate(time_s=float(result.time[i]), snr=float(result.snr[i]),
                  chi2_reduced=shape_veto_chi2(result, float(result.time[i])))
        for i in peak_idx
    ]

    if max_chi2_reduced is not None:
        candidates = [c for c in candidates if c.chi2_reduced <= max_chi2_reduced]

    candidates.sort(key=lambda c: c.snr, reverse=True)
    return candidates
