# kbo_occultation/filtering.py
"""
High-pass filtering of the fast light curve, to remove the correlated,
sub-few-Hz noise excess found in inspect_real_noise_2025_12_16.py without
touching the frequencies above it, where the occultation signal itself
lives.

Deliberately independent of dc_combine.py's despiking/detrending: this only
needs a uniformly-sampled light curve and a cutoff frequency, so it can be
applied instead of, before, or after any DC-based correction, and doesn't
need to change if a better noise-removal method replaces despiking later.

The one thing that *does* couple the two: despike_lightcurve_with_dc's
default (fill="remove") leaves gaps in the time grid, and both filters here
assume uniform sampling like any FFT/IIR method would. If combining the
two, high-pass filter first (on the still-gapless raw or DC-corrected
signal) and despike afterward.
"""

import numpy as np
from scipy.signal import butter, sosfiltfilt
from scipy.ndimage import median_filter


def highpass_fft(signal: np.ndarray, fs: float, cutoff_hz: float) -> np.ndarray:
    """
    Brick-wall high-pass filter: zero every rfft bin below cutoff_hz (and
    the signal's mean, which is bin 0) and invert, then add the mean back
    so the output stays in the original signal's units/offset rather than
    oscillating around zero.

    Exact in the frequency domain -- no rolloff/order to tune -- but a hard
    cutoff can ring in the time domain around sharp features (a bin just
    above cutoff_hz is untouched, one just below is zeroed outright). See
    highpass_butterworth for a smoother alternative.
    """
    signal = np.asarray(signal)
    mean = np.mean(signal)

    freq_hz = np.fft.rfftfreq(len(signal), d=1.0 / fs)
    spectrum = np.fft.rfft(signal - mean)
    spectrum[freq_hz < cutoff_hz] = 0.0

    return np.fft.irfft(spectrum, n=len(signal)) + mean


def highpass_butterworth(signal: np.ndarray, fs: float, cutoff_hz: float, order: int = 4) -> np.ndarray:
    """
    Zero-phase Butterworth high-pass (sosfiltfilt applies the filter
    forward and backward, so there's no time shift and no ringing
    precursor sharper than the filter's own smooth rolloff) -- a gentler
    alternative to highpass_fft's brick-wall cutoff, at the cost of a
    transition band around cutoff_hz rather than an exact one.

    As in highpass_fft, the mean is subtracted before filtering (a
    high-pass response is 0 at DC anyway, but this avoids edge-padding
    transients tied to a large fixed offset) and added back after.
    """
    signal = np.asarray(signal)
    mean = np.mean(signal)

    sos = butter(order, cutoff_hz, btype="highpass", fs=fs, output="sos")
    return sosfiltfilt(sos, signal - mean) + mean


def highpass_lightcurve(lc, cutoff_hz: float, method: str = "butterworth", **kwargs):
    """
    Apply a high-pass filter to a LightCurve's signal.

    Parameters
    ----------
    lc : LightCurve
        Must be uniformly sampled -- true of every LightCurve.from_stat_binary*
        loader in this package and of dc_combine's fft_swap/detrend outputs,
        but *not* despike_lightcurve_with_dc's default fill="remove" output
        (see module docstring for the ordering that avoids this).
    cutoff_hz : float
        High-pass cutoff, in Hz.
    method : {"butterworth", "fft"}
        See highpass_butterworth / highpass_fft.
    **kwargs
        Passed through to the chosen method (e.g. order= for butterworth).
    """
    from .photometry import LightCurve  # local import: photometry imports this module too

    fs = 1.0 / lc.dt
    if method == "butterworth":
        filtered = highpass_butterworth(lc.signal, fs, cutoff_hz, **kwargs)
    elif method == "fft":
        filtered = highpass_fft(lc.signal, fs, cutoff_hz, **kwargs)
    else:
        raise ValueError(f"method must be 'butterworth' or 'fft', got {method!r}")

    meta = dict(lc.meta)
    meta.update({"highpass_cutoff_hz": cutoff_hz, "highpass_method": method})
    return LightCurve(lc.time, filtered, meta=meta)


def remove_outliers_lightcurve(lc, threshold_sigma: float = 9.0, window_samples: int = 21,
                               polarity: str = "positive"):
    """
    Remove isolated single-sample spikes from a LightCurve's signal.

    Self-contained (no DC report needed, unlike dc_combine.despike_lightcurve_with_dc).
    Targets the narrow upward spikes (typically one sample) seen in the fast
    curves -- cosmic hits and readout glitches -- *without* touching real
    structure: the broad wind/tracking dips or the ms-scale occultation dips the
    search pipeline is meant to find.

    Two choices make it specific to those spikes rather than a general sigma clip:

    - The baseline is a *short* rolling median (``window_samples``, a handful of
      samples), so it hugs the real signal and only a point that sticks out from
      its immediate neighbours has a large residual. A wide window would instead
      straddle broad features and flag their extrema.
    - Only ``polarity="positive"`` residuals are flagged by default, so downward
      dips (real occultations, wind valleys) are left completely alone.

    Flagged samples are *replaced by the local baseline value* (the neighbours'
    median), which keeps the array length and uniform time grid intact -- so the
    output stays valid input to highpass_lightcurve, PSD estimation and the
    matched filter, unlike dropping samples.

    Parameters
    ----------
    lc : LightCurve
        Must be uniformly sampled (true of every LightCurve.from_stat_binary*
        loader).
    threshold_sigma : float
        Residuals beyond this many robust (MAD-based) sigmas are flagged.
    window_samples : int
        Width of the short rolling-median baseline, in samples (forced odd,
        minimum 3). Keep it small -- a handful of samples -- so only isolated
        spikes, not broad features, are caught.
    polarity : {"positive", "negative", "both"}
        Which side to clip. "positive" (default) flags only upward spikes;
        "negative" only downward ones; "both" is a symmetric clip.

    Returns
    -------
    LightCurve
        New LightCurve on the same time grid with flagged samples replaced by
        the local baseline. ``meta`` gains ``outlier_threshold_sigma``,
        ``outlier_window_samples``, ``outlier_polarity`` and ``outlier_fraction``.
    """
    from .photometry import LightCurve  # local import: photometry imports this module too
    from .matched_filter import robust_sigma

    if polarity not in ("positive", "negative", "both"):
        raise ValueError(f"polarity must be 'positive', 'negative' or 'both', got {polarity!r}")

    signal = np.asarray(lc.signal)

    window_samples = max(3, int(window_samples))
    if window_samples % 2 == 0:
        window_samples += 1

    baseline = median_filter(signal, size=window_samples, mode="nearest")
    residual = signal - baseline
    sigma = robust_sigma(residual)
    if sigma == 0:
        # If the residual is mostly exactly zero (heavily quantized or very
        # stable stretches), the MAD-based robust sigma degenerates to zero and
        # any nonzero residual would trivially "exceed" it. Fall back to the
        # plain std, which stays well-defined as long as there is real scatter.
        # (Same guard as dc_combine.flag_dc_excursions.)
        sigma = float(np.std(residual))

    if polarity == "positive":
        bad = residual > threshold_sigma * sigma
    elif polarity == "negative":
        bad = residual < -threshold_sigma * sigma
    else:
        bad = np.abs(residual) > threshold_sigma * sigma

    cleaned = signal.copy()
    cleaned[bad] = baseline[bad]

    meta = dict(lc.meta)
    meta.update({
        "outlier_threshold_sigma": threshold_sigma,
        "outlier_window_samples": window_samples,
        "outlier_polarity": polarity,
        "outlier_fraction": float(np.mean(bad)),
    })
    return LightCurve(lc.time, cleaned, meta=meta)
