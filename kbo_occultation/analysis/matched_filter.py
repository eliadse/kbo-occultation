# kbo_occultation/analysis/matched_filter.py
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
from typing import List, Literal, Optional, Tuple

import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.signal import correlate, find_peaks

from ..stats import robust_sigma  # re-exported: robust_sigma used to live here


@dataclass
class MatchedFilterResult:
    time: np.ndarray
    snr: np.ndarray
    # Local (per-sample) robust noise level, aligned with data_time /
    # residual_data -- see compute_baseline_and_sigma.
    sigma: np.ndarray
    baseline: np.ndarray
    dt_s: float
    template_duration_s: float
    # Kept around so shape_veto_chi2 can re-fit the template against any
    # candidate's native-resolution data window without redoing the
    # detrending/resampling step.
    data_time: np.ndarray
    residual_data: np.ndarray
    template_delta: np.ndarray
    # Per-lag noise of the filter output in flux units, i.e. the white
    # per-sample sigma that would give the same SNR (snr = corr /
    # template_norm / amp_sigma), aligned with `time`. For
    # normalization="empirical" this is measured on the output itself;
    # otherwise it is the white-noise value sigma_center (see
    # sliding_matched_filter_snr).
    amp_sigma: Optional[np.ndarray] = None
    template_norm: float = 1.0
    normalization: str = "empirical"
    # empirical_filter_sigma method behind amp_sigma (normalization="empirical" only)
    sigma_method: Optional[str] = None


@dataclass
class Candidate:
    time_s: float
    snr: float
    chi2_reduced: float


BaselineMethod = Literal["median", "clipped_triangle"]


def _clipped_triangle(x: np.ndarray, window: int, clip_k: float = 3.0) -> np.ndarray:
    """
    Robust but linear smoothing of ``x`` over ``window`` samples: points more
    than ``clip_k`` robust sigmas from the rolling median are replaced by that
    median, then a triangular kernel of full width ``window`` (two running
    means of window/2) is applied.
    """
    med = median_filter(x, size=window, mode="reflect")
    dev = x - med
    spread = 1.4826 * median_filter(np.abs(dev), size=window, mode="reflect")
    kept = np.where(np.abs(dev) > clip_k * spread, med, x)
    half = max(window // 2, 1)
    return uniform_filter1d(uniform_filter1d(kept, size=half, mode="reflect"), size=half, mode="reflect")


def compute_baseline_and_sigma(time: np.ndarray, signal: np.ndarray, template_duration_s: float,
                                detrend_window_s: Optional[float] = None,
                                decimate="auto",
                                method: BaselineMethod = "median") -> Tuple[np.ndarray, np.ndarray]:
    """
    Rolling baseline and a local,
    per-sample robust noise level of the resulting residual, factored
    out of sliding_matched_filter_snr so it can be computed once and reused
    (e.g. injection_batch.py caches this per noise-reduction method
    instead of recomputing it for every trial).

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
    decimate : "auto", int
        Average blocks of `decimate` samples, roll the median over the
        block means, and linearly interpolate back to full resolution.
        The baseline is smooth on scales far below the window by
        construction, and block-averaging first keeps the decimated
        median's statistical scatter at the exact filter's level (a
        plain strided median over window/stride points would be
        sqrt(stride) noisier), so this matches the exact result to
        ~1e-3 relative on %-level noise while cutting the cost by
        ~stride^2 (the exact filter is O(N * window)). A short dip or
        spike only contaminates a few block means, which the median
        across the window still rejects. "auto" picks
        stride = window // 64; pass an int for an explicit stride, or
        1 for the exact (slow) full-resolution median.
    method : {"median", "clipped_triangle"}
        ``"median"``: rolling median (scipy.ndimage.median_filter). Robust,
        but on red noise it moves in steps (it switches between values), and
        the template, whose mean is not zero, responds to steps: the filter
        output of Gaussian noise gets heavy tails (excess kurtosis ~0.5-1).
        ``"clipped_triangle"``: values more than 3 robust sigmas from the
        rolling median (a dip, a glitch, a DC excursion) are replaced by that
        median, then a triangular kernel of full width the window smooths the
        rest. Linear wherever nothing is clipped, so Gaussian noise stays
        Gaussian; the triangle's weak frequency sidelobes let less slow wander
        into the filter output than a plain running mean
        (notebooks/baseline_study.ipynb).

    Returns
    -------
    baseline : ndarray, same shape as signal.
    sigma : ndarray, same shape as signal. Local robust noise level
        (MAD-based) of (signal/baseline - 1), so a night whose noise
        drifts (elevation, transparency, sky brightness) gets a sigma
        that tracks it per baseline block instead of one number for the
        whole series. Estimated per decimated block and interpolated
        back to full resolution, mirroring the baseline.

    Caveat -- this sigma is a *white-noise* level. It is the MAD of the
    scatter within blocks (centred on each block's own median), so it sees
    only noise faster than ~1 / block duration. On red noise the decimated
    path ("auto") reads ~30 % lower than the exact one (decimate=1), and
    either one understates the noise in the template's band by a large
    factor (SNR sigma 3.6-5.3 instead of 1 in the test surrogates). The
    honest normalisation is sliding_matched_filter_snr(normalization=
    "empirical"), which is insensitive to this choice.
    """
    time = np.asarray(time)
    signal = np.asarray(signal)
    dt = float(np.median(np.diff(time)))

    if detrend_window_s is None:
        detrend_window_s = 20.0 * template_duration_s
    M = int(round(template_duration_s / dt))
    window_samples = max(int(round(detrend_window_s / dt)), M)

    if decimate == "auto":
        stride = max(1, window_samples // 64)
    else:
        stride = max(1, int(decimate))

    if method not in ("median", "clipped_triangle"):
        raise ValueError(f"Unknown baseline method {method!r}")

    if stride == 1:
        if method == "median":
            baseline = median_filter(signal, size=window_samples, mode="reflect")
        else:
            baseline = _clipped_triangle(signal, window_samples)
        residual = signal / baseline - 1.0
        # After detrending the residual is ~zero-median locally, so a
        # rolling median of |residual| over the window is the local MAD.
        sigma = 1.4826 * median_filter(np.abs(residual), size=window_samples,
                                       mode="reflect")
    else:
        n_blocks = len(signal) // stride
        sig_blocks = signal[:n_blocks * stride].reshape(-1, stride).mean(axis=1)
        t_blocks = time[:n_blocks * stride].reshape(-1, stride).mean(axis=1)
        window_dec = max(window_samples // stride, 3)
        if method == "median":
            baseline_dec = median_filter(sig_blocks, size=window_dec, mode="reflect")
        else:
            baseline_dec = _clipped_triangle(sig_blocks, window_dec)
        baseline = np.interp(time, t_blocks, baseline_dec)

        residual = signal / baseline - 1.0
        # Per-block MAD, then median-smoothed over the window (rejecting
        # blocks contaminated by a dip/spike, same argument as the
        # baseline) and interpolated back to full resolution.
        res_blocks = residual[:n_blocks * stride].reshape(-1, stride)
        block_mad = np.median(
            np.abs(res_blocks - np.median(res_blocks, axis=1, keepdims=True)),
            axis=1,
        )
        sigma_dec = 1.4826 * median_filter(block_mad, size=window_dec, mode="reflect")
        sigma = np.interp(time, t_blocks, sigma_dec)

    return baseline, sigma


Normalization = Literal["local_block", "global", "empirical"]
SigmaMethod = Literal["sliding", "blocks"]


def empirical_filter_sigma(amplitude: np.ndarray, window_samples: int,
                           block_samples: Optional[int] = None,
                           method: SigmaMethod = "blocks") -> np.ndarray:
    """
    Rolling robust (1.4826 * median|x|) scatter of the matched-filter output
    itself, per lag. Because the output of a filter on red noise is strongly
    correlated, this -- unlike sigma_data / sqrt(sum T^2) -- includes the
    noise power inside the template's own band.

    ``method="sliding"``: the median of |x| over a window of
    ``window_samples`` centred on every lag (evaluated every
    ``window_samples // 320`` lags and interpolated: the output is smooth on
    the template scale, so this loses nothing). No block edges, so the
    estimate does not depend on where the series starts.

    ``method="blocks"``: the output is measured block by block (blocks of
    ``block_samples``, default ``window_samples // 8``), and the per-block
    values are median-smoothed over ``window_samples`` and interpolated back
    to every lag. On red noise the filter output is correlated over about
    one block, so each block holds ~1 independent value and the estimate
    depends on the block edges by tens of percent, and the median of the
    skewed block values reads ~7 % low (notebooks/sigma_estimator_study.ipynb).

    An event occupies about one template length, a small part of the window,
    and the median (50 % breakdown) rejects it, which is the "exclude +-1
    template around the peak" of the plan done implicitly. With "sliding" a
    strong event still raises the estimate at its own lag by ~1.3x the
    fraction of the window it covers (~4-8 % for 1-5 %); "blocks" confines it
    to one block value.
    """
    if method == "sliding":
        return _sliding_filter_sigma(np.asarray(amplitude), window_samples)
    if method != "blocks":
        raise ValueError(f"Unknown sigma method {method!r}")
    n = len(amplitude)
    block = max(int(block_samples or max(window_samples // 8, 8)), 8)
    n_blocks = n // block
    if n_blocks < 3:
        return np.full(n, 1.4826 * np.median(np.abs(amplitude)))
    chunks = amplitude[:n_blocks * block].reshape(n_blocks, block)
    # Median of |amp| about zero, NOT about each block's own median: the
    # output is already detrended, and re-centering per block would throw
    # away the slow wander of correlated noise, which is exactly the
    # excess this estimator exists to capture.
    mad = np.median(np.abs(chunks), axis=1)
    win_blocks = max(window_samples // block, 3) | 1
    sigma_blocks = 1.4826 * median_filter(mad, size=win_blocks, mode="reflect")
    centers = (np.arange(n_blocks) + 0.5) * block
    return np.interp(np.arange(n), centers, sigma_blocks)


def _sliding_filter_sigma(amplitude: np.ndarray, window_samples: int) -> np.ndarray:
    """empirical_filter_sigma(method="sliding"): centred rolling median of |x|."""
    n = len(amplitude)
    step = max(1, int(window_samples) // 320)
    sub = np.abs(amplitude[::step])
    if len(sub) < 3:
        return np.full(n, 1.4826 * np.median(np.abs(amplitude)))
    size = max(int(window_samples) // step, 3) | 1
    med = median_filter(sub, size=size, mode="reflect")
    return np.interp(np.arange(n), np.arange(len(sub)) * step, 1.4826 * med)


def sliding_matched_filter_snr(time: np.ndarray, signal: np.ndarray,
                                template_time: np.ndarray, template_intensity: np.ndarray,
                                detrend_window_s: Optional[float] = None,
                                baseline: Optional[np.ndarray] = None,
                                sigma: Optional[np.ndarray] = None,
                                normalization: Normalization = "empirical",
                                amp_sigma: Optional[np.ndarray] = None,
                                sigma_window_s: Optional[float] = None,
                                sigma_method: SigmaMethod = "blocks",
                                baseline_method: BaselineMethod = "median") -> MatchedFilterResult:
    """
    Slide (template_time, template_intensity) a noiseless simulated diffraction pattern 
    in time units, e.g. via detectability.spatial_to_time, across (time, signal) 
    and compute the matched-filter SNR at every lag.

    Parameters
    ----------
    time, signal : array_like
        Light curve to search (real or injected-real), on a uniform time grid.
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
    normalization : {"empirical", "local_block", "global"}
        How the correlation is turned into an SNR.

        - ``"empirical"`` (default): divide by the rolling robust scatter
          of the filter output itself (empirical_filter_sigma). Honest for
          correlated noise: SNR on pure noise has sigma ~ 1 whatever the
          noise spectrum.
        - ``"local_block"``: divide by template_norm * the per-block MAD
          of the residual (``sigma``). This assumes white noise and sees
          only noise above ~1/(block duration) -- optimistic on red
          noise. Kept for regression.
        - ``"global"``: as local_block with one robust sigma for the
          whole series.
    sigma_window_s : float, optional
        Window over which ``normalization="empirical"`` measures the scatter
        of the filter output. Defaults to the detrend window (20x the
        template). On noise coherent over many seconds a window this short
        holds only a few independent samples and the sigma itself fluctuates
        by tens of percent along the run, which fattens the SNR tails; use a
        longer window (a few minutes, or the whole run) when the noise is
        stationary enough.
    sigma_method : {"blocks", "sliding"}
        ``empirical_filter_sigma`` method for ``normalization="empirical"``.
        "sliding" (a centred rolling median, no block edges) is the less
        noisy and unbiased one; "blocks" is the original.
    baseline_method : {"median", "clipped_triangle"}
        ``compute_baseline_and_sigma`` method, when the baseline is computed
        here (i.e. not passed in).
    amp_sigma : ndarray, optional
        Precomputed per-lag filter-output sigma (``MatchedFilterResult.
        amp_sigma`` of a search on the same template), length
        ``len(signal) - M + 1``. Overrides ``normalization``. Used by the
        injection Monte Carlo to normalise every short trial slice with the
        sigma measured once on the whole uninjected series: measured on the
        slice itself the estimate scatters ~10 % trial to trial.

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
        baseline, sigma = compute_baseline_and_sigma(time, signal, template_duration_s, detrend_window_s,
                                                     method=baseline_method)
    norm_signal = signal / baseline

    # --- matched filter: FFT cross-correlation of the detrended data with the template variation ---
    residual_data = norm_signal - 1.0
    template_delta = template_resampled - 1.0
    template_norm = float(np.sqrt(np.sum(template_delta**2)))

    correlation = correlate(residual_data, template_delta, mode="valid", method="fft")

    # correlation[i] is the match for the window signal[i:i+M]; tag it
    # at that window's center time. sigma varies slowly (block/window
    # scale >> template), so the local noise at each window's center is
    # the right normaliser for that lag.
    center_offset = M // 2
    n_lag = len(correlation)
    if amp_sigma is not None:
        amp_sigma = np.asarray(amp_sigma, dtype=float)
        if len(amp_sigma) != n_lag:
            raise ValueError(f"amp_sigma has {len(amp_sigma)} lags, expected {n_lag}")
        normalization = "supplied"
    elif normalization == "empirical":
        window_s = sigma_window_s or detrend_window_s or 20.0 * template_duration_s
        window_samples = max(int(round(window_s / dt)), M)
        amp_sigma = empirical_filter_sigma(correlation / template_norm, window_samples, method=sigma_method)
    elif normalization == "global":
        amp_sigma = np.full(n_lag, robust_sigma(residual_data))
    elif normalization == "local_block":
        amp_sigma = sigma[center_offset:center_offset + n_lag]
    else:
        raise ValueError(f"Unknown normalization {normalization!r}")
    snr = correlation / (template_norm * amp_sigma)

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
        template_delta=template_delta,
        amp_sigma=amp_sigma,
        template_norm=template_norm,
        normalization=normalization,
        sigma_method=sigma_method if normalization == "empirical" else None,
    )


def shape_veto_chi2(result: MatchedFilterResult, candidate_time_s: float,
                    weighting: Literal["local_block", "empirical"] = "local_block") -> float:
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
    weighting : {"local_block", "empirical"}
        Per-sample noise used to weight the residuals. ``"local_block"``
        uses result.sigma (white-noise MAD). ``"empirical"`` uses the
        effective per-sample sigma that reproduces the measured scatter
        of the filter output, ``amp_sigma`` (flux units; the sigma a white
        series would need to give the same SNR): for white noise this equals sigma, for red noise it is larger, so chi2/dof
        is no longer inflated by correlated noise that the SNR already
        accounts for. Use calibrate_chi2_threshold to set the cut.

    Returns
    -------
    float
        chi2 / dof, dof = (number of samples in the template window - 1).
    """
    M = len(result.template_delta)
    center_idx = int(np.searchsorted(result.data_time, candidate_time_s))
    start = max(0, min(center_idx - M // 2, len(result.residual_data) - M))
    window = result.residual_data[start:start + M]

    template_norm_sq = float(np.sum(result.template_delta**2))
    amplitude = float(np.sum(window * result.template_delta) / template_norm_sq)

    residuals = window - amplitude * result.template_delta
    if weighting == "empirical":
        if result.amp_sigma is None:
            raise ValueError("result carries no amp_sigma; recompute with sliding_matched_filter_snr")
        lag = int(np.clip(center_idx - M // 2, 0, len(result.amp_sigma) - 1))
        local_sigma = result.amp_sigma[lag]
    elif weighting == "local_block":
        local_sigma = result.sigma[start:start + M]
    else:
        raise ValueError(f"Unknown weighting {weighting!r}")
    dof = M - 1
    return float(np.sum((residuals / local_sigma)**2) / dof)


def calibrate_chi2_threshold(result: MatchedFilterResult, quantile: float = 0.99,
                             weighting: Literal["local_block", "empirical"] = "empirical",
                             n_samples: int = 400, exclude_snr: float = 3.0,
                             rng: Optional[np.random.Generator] = None) -> Tuple[float, float]:
    """
    Calibrate max_chi2_reduced on off-source data: evaluate the shape veto
    at random lags whose SNR is below ``exclude_snr`` (so real events are
    not part of the calibration sample) and return
    (median chi2/dof, quantile chi2/dof). The median should be ~1 on pure
    noise with the right weighting; the quantile is a data-driven
    ``max_chi2_reduced``.
    """
    rng = rng or np.random.default_rng(0)
    ok = np.flatnonzero(np.abs(result.snr) < exclude_snr)
    if ok.size == 0:
        raise ValueError("no off-source lags below exclude_snr")
    picks = rng.choice(ok, size=min(n_samples, ok.size), replace=False)
    chi = np.array([shape_veto_chi2(result, float(result.time[i]), weighting) for i in picks])
    return float(np.median(chi)), float(np.quantile(chi, quantile))


def find_candidates(result: MatchedFilterResult, snr_threshold: float = 5.0,
                     min_separation_s: Optional[float] = None,
                     max_chi2_reduced: Optional[float] = None,
                     chi2_weighting: Literal["local_block", "empirical"] = "local_block"
                     ) -> List[Candidate]:
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
                  chi2_reduced=shape_veto_chi2(result, float(result.time[i]), chi2_weighting))
        for i in peak_idx
    ]

    if max_chi2_reduced is not None:
        candidates = [c for c in candidates if c.chi2_reduced <= max_chi2_reduced]

    candidates.sort(key=lambda c: c.snr, reverse=True)
    return candidates
