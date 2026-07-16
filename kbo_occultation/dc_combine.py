# kbo_occultation/dc_combine.py
"""
Combine the fast, high-cadence correlator signal (std_dev^2, ~1.9 kHz) with
the independently-recorded slow DC report (~1 Hz) into a single composite
light curve.

Motivation (see examples/inspect_real_noise_2025_12_16.py): the fast
signal's noise is dominated by a correlated, sub-few-Hz component that
block-averaging cannot remove, and that component is highly coherent with
the slow DC channel -- i.e. the same real (tracking/wind/transparency)
systematic is visible, less noisily, in the slow monitoring. Rather than
discard the fast signal's high-frequency content (where an occultation dip
actually lives), we swap in the DC report's own frequency content below a
cutoff and keep the fast signal's above it.

Direct port of the `combined_dc_reports_with_std_dev_2` prototype developed
in kbo_occultation_rates' notebooks/check_fast_photometry_data.ipynb, split
into a generic array-level function (combine_fast_and_dc) and a
LightCurve/DC-report-file convenience (combine_lightcurve_with_dc).
"""

import pickle
from typing import Optional, Tuple

import numpy as np
from scipy.ndimage import median_filter
from astropy.time import Time

from .photometry import LightCurve
from .matched_filter import robust_sigma


def load_dc_series(dc_report_path: str, telescope: str = "m2", pixel: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load a slow DC report pickle and return one pixel's series as
    (time_unix_s, value).

    Parameters
    ----------
    dc_report_path : str
        Path to a dc_report.pkl, as produced by magic_spysii's
        DcReportReader (e.g. data/observations/DCs/2025_12_16/dc_report.pkl).
        Unpickling requires the `magic_spysii` package to be importable.
    telescope : str
        Attribute name on the DcReportReader holding this telescope's
        report, e.g. "m2" (MAGIC-II), "m1" (MAGIC-I), "lst1" (LST1).
    pixel : int
        Column of the `dcs` array to read -- 0 is the on-target pixel,
        1/2 are background pixels (see draw_dcs()-style plots in the
        original notebook).
    """
    with open(dc_report_path, "rb") as f:
        dc_report = pickle.load(f)

    tel_data = getattr(dc_report, telescope)
    if len(tel_data["times"]) == 0:
        raise ValueError(f"No DC data for telescope {telescope!r} in {dc_report_path}")

    dc_time_s = Time(tel_data["times"], format="isot", scale="utc").unix
    dc_value = np.asarray(tel_data["dcs"][:, pixel])
    return dc_time_s, dc_value


def interpolate_dc_to_fast(fast_time_s: np.ndarray, dc_time_s: np.ndarray, dc_value: np.ndarray) -> np.ndarray:
    """
    Linearly interpolate a slow DC series onto the fast light curve's time
    grid (both must be unix seconds, or any other shared linear time base).

    Raises if `fast_time_s` isn't covered by `dc_time_s` to more than one DC
    sampling interval on either end -- outside that, np.interp would
    silently clamp to the first/last DC value rather than extrapolate,
    which is fine right at the edge but a sign of a real time-base mismatch
    if it stretches further than that.
    """
    fast_time_s = np.asarray(fast_time_s)
    dc_dt_s = float(np.median(np.diff(dc_time_s)))

    if fast_time_s.min() < dc_time_s.min() - dc_dt_s or fast_time_s.max() > dc_time_s.max() + dc_dt_s:
        raise ValueError(
            "Fast light curve extends beyond the DC report's time coverage by more than one "
            "DC sampling interval -- check that both come from the same observation."
        )

    return np.interp(fast_time_s, dc_time_s, dc_value)


def combine_fast_and_dc(time_s: np.ndarray, fast_signal: np.ndarray, dc_on_fast: np.ndarray,
                         frequency_cut_hz: float, return_spectra: bool = False):
    """
    Frequency-domain hybrid of a fast signal and a slow DC series already
    resampled onto the same uniform time grid: below frequency_cut_hz, use
    the DC report's own spectral content; at and above it, use the fast
    signal's. Both are mean-normalised before combining, and the result is
    re-scaled back to the DC's absolute level.

    Caution -- this is a hard, bin-wise spectral swap (equivalent to
    filtering.highpass_fft's brick-wall cutoff applied to the fast signal,
    summed with an equally brick-wall low-pass of the DC series), not a
    smooth rolloff like filtering.highpass_butterworth. As that module's
    own docstring notes, a brick-wall cutoff can ring in the time domain
    around sharp features -- and an injection-recovery test
    (kbo_occultation.injection_batch.run_injection_monte_carlo, see
    examples/injection_montecarlo_2025_12_16.py) confirmed this matters in
    practice: for a real KBO occultation template whose own timescale
    (~0.3s here) is only ~3x this function's default cutoff period (1s),
    the combined signal's fine structure can diverge enough from the raw
    template's that a shape-consistency check
    (matched_filter.shape_veto_chi2) rejects a genuine, high-SNR event
    (measured directly: SNR~165 but chi2/dof~18 for one R=300m trial).
    Rebuilding the comparison template by running this same function on a
    flat reference (matched_filter/injection_batch's approach for a true
    LTI filter like highpass) only partly compensates (chi2/dof dropped to
    ~5, still short of passing) with a few seconds of context on each side
    -- consistent with the ringing extending further than that, unlike the
    fast-decaying, near-local response of a Butterworth filter. Net effect:
    this correction can distort real event shape for occultations whose
    duration approaches 1/frequency_cut_hz, in a way that's harder to
    correct for in a matched-filter search than highpass's cleaner,
    faster-decaying distortion is.

    Parameters
    ----------
    time_s : array_like
        Uniform time grid shared by `fast_signal` and `dc_on_fast`.
    fast_signal : array_like
        High-cadence signal (e.g. LightCurve.signal, std_dev^2).
    dc_on_fast : array_like
        Slow DC series already interpolated onto `time_s` (see
        interpolate_dc_to_fast).
    frequency_cut_hz : float
        Crossover frequency: rfft bins with frequency < frequency_cut_hz
        come from `dc_on_fast`; bins at or above it come from `fast_signal`.
        The PSD/coherence checks in inspect_real_noise_2025_12_16.py put the
        correlated noise excess below ~1-2 Hz, so that's a reasonable
        starting point -- but it should stay well below the frequency of
        the occultation signal you're trying to preserve.
    return_spectra : bool
        If True, also return the individual and combined rfft spectra
        (mean-normalised), for diagnosing the choice of frequency_cut_hz.

    Returns
    -------
    combined : ndarray
        Same length as the inputs, in the DC's absolute units.
    spectra : dict, optional
        {'freq_hz', 'fft_dc', 'fft_fast', 'fft_combined'} -- only if
        return_spectra.
    """
    time_s = np.asarray(time_s)
    fast_signal = np.asarray(fast_signal)
    dc_on_fast = np.asarray(dc_on_fast)

    dt_s = float(np.median(np.diff(time_s)))
    freq_hz = np.fft.rfftfreq(len(time_s), d=dt_s)

    mean_dc = np.mean(dc_on_fast)
    mean_fast = np.mean(fast_signal)

    fft_dc = np.fft.rfft(dc_on_fast / mean_dc)
    fft_fast = np.fft.rfft(fast_signal / mean_fast)

    fft_combined = fft_dc.copy()
    above_cut = freq_hz >= frequency_cut_hz
    fft_combined[above_cut] = fft_fast[above_cut]

    combined = np.fft.irfft(fft_combined, n=len(time_s)) * mean_dc

    if return_spectra:
        return combined, {
            "freq_hz": freq_hz,
            "fft_dc": fft_dc,
            "fft_fast": fft_fast,
            "fft_combined": fft_combined,
        }
    return combined


def combine_lightcurve_with_dc(lc: LightCurve, dc_report_path: str, frequency_cut_hz: float,
                                telescope: str = "m2", pixel: int = 0,
                                return_spectra: bool = False) -> LightCurve:
    """
    End-to-end convenience: load a DC report, interpolate it onto `lc`'s
    time grid, and return a new LightCurve whose signal is the frequency-
    domain hybrid from combine_fast_and_dc -- see that function's
    docstring for a real, measured caveat about brick-wall ringing
    distorting occultation shape for a matched-filter search.

    Parameters
    ----------
    lc : LightCurve
        Fast light curve, e.g. from LightCurve.from_stat_binary_all(...)[channel].
    dc_report_path, telescope, pixel : see load_dc_series.
    frequency_cut_hz, return_spectra : see combine_fast_and_dc.
    """
    dc_time_s, dc_value = load_dc_series(dc_report_path, telescope=telescope, pixel=pixel)
    dc_on_fast = interpolate_dc_to_fast(lc.time, dc_time_s, dc_value)

    result = combine_fast_and_dc(lc.time, lc.signal, dc_on_fast, frequency_cut_hz, return_spectra=return_spectra)
    combined, spectra = result if return_spectra else (result, None)

    meta = dict(lc.meta)
    meta.update({
        "dc_report_path": dc_report_path,
        "dc_telescope": telescope,
        "dc_pixel": pixel,
        "dc_frequency_cut_hz": frequency_cut_hz,
    })
    combined_lc = LightCurve(lc.time, combined, meta=meta)

    if return_spectra:
        return combined_lc, spectra
    return combined_lc


# ─── Alternative uses of the DC: a smoothed trend, and a bad-data mask ────
#
# combine_fast_and_dc (above) imports the DC report's *raw* per-sample
# values below the cutoff -- which turned out to make 2025-12-16 worse, not
# better, because that channel's own sample-to-sample scatter (~4-8%
# fractional, at ~1 Hz) is larger than the fast signal's already-elevated
# low-frequency scatter. The two functions below instead treat the DC as
# (a) a noisy proxy for a *slow trend* worth smoothing before use, and (b)
# a source of a bad-data mask for discrete dropouts (wind gusts, tracking
# glitches), rather than a frequency-domain amplitude reference.

def smooth_dc_series(dc_time_s: np.ndarray, dc_value: np.ndarray, window_s: float) -> np.ndarray:
    """
    Rolling-median smoothing of a DC series at its own (native, ~1 Hz)
    cadence, to average out this channel's own sample-to-sample scatter
    while preserving genuine multi-second trends (tracking/wind/
    transparency) -- the only part of the DC series that should be
    imported into the fast signal (see detrend_lightcurve_with_dc).
    """
    dc_dt_s = float(np.median(np.diff(dc_time_s)))
    window_samples = max(3, int(round(window_s / dc_dt_s)))
    if window_samples % 2 == 0:
        window_samples += 1
    return median_filter(dc_value, size=window_samples, mode="nearest")


def detrend_lightcurve_with_dc(lc: LightCurve, dc_report_path: str, telescope: str = "m2", pixel: int = 0,
                                smooth_window_s: float = 10.0) -> LightCurve:
    """
    Divide the fast light curve by a smoothed, mean-normalised DC trend
    (smooth_dc_series), removing the shared slow systematic without
    importing the DC channel's own sample-to-sample noise -- unlike
    combine_lightcurve_with_dc, which uses the DC's raw values.

    Parameters
    ----------
    lc : LightCurve
        Fast light curve.
    dc_report_path, telescope, pixel : see load_dc_series.
    smooth_window_s : float
        Median-filter window (seconds) applied to the DC series before
        it's used as a trend. Should be long enough to average out the
        DC's own per-sample scatter but short enough to still track real
        transparency/tracking drift -- 10s is a starting point, well above
        the DC's own ~1s cadence and well below a 15-minute run.
    """
    dc_time_s, dc_value = load_dc_series(dc_report_path, telescope=telescope, pixel=pixel)

    # Restrict to this light curve's own observing window before smoothing
    # (matches despike_lightcurve_with_dc): the DC report runs continuously
    # through the whole night, across every target and the idle/slewing
    # gaps between them, so smoothing the full series is both wasteful and
    # (for a long enough window relative to a gap) a route for a
    # neighbouring target's completely different DC level to leak into this
    # one's edge samples.
    in_window = (dc_time_s >= lc.time.min() - smooth_window_s) & (dc_time_s <= lc.time.max() + smooth_window_s)
    dc_time_s, dc_value = dc_time_s[in_window], dc_value[in_window]

    dc_smooth = smooth_dc_series(dc_time_s, dc_value, smooth_window_s)
    trend_on_fast = interpolate_dc_to_fast(lc.time, dc_time_s, dc_smooth)
    trend_norm = trend_on_fast / np.mean(trend_on_fast)

    meta = dict(lc.meta)
    meta.update({
        "dc_report_path": dc_report_path,
        "dc_telescope": telescope,
        "dc_pixel": pixel,
        "dc_smooth_window_s": smooth_window_s,
        "correction": "smoothed_dc_detrend",
    })
    return LightCurve(lc.time, lc.signal / trend_norm, meta=meta)


def flag_dc_excursions(dc_time_s: np.ndarray, dc_value: np.ndarray, baseline_window_s: float = 60.0,
                        threshold_sigma: float = 4.0) -> np.ndarray:
    """
    Flag DC samples that deviate from a rolling-median baseline by more
    than threshold_sigma robust sigmas -- the discrete wind/tracking dips
    the original notebook's half-built remove_valleys_dc was chasing --
    rather than trying to correct their amplitude.

    Returns a boolean mask on dc_time_s/dc_value, True where flagged.
    """
    dc_dt_s = float(np.median(np.diff(dc_time_s)))
    window_samples = max(3, int(round(baseline_window_s / dc_dt_s)))
    if window_samples % 2 == 0:
        window_samples += 1
    baseline = median_filter(dc_value, size=window_samples, mode="nearest")
    residual = dc_value - baseline
    sigma = robust_sigma(residual)
    if sigma == 0:
        # The DC readout is coarsely quantized (~0.01 steps) and often
        # changes slowly relative to that step, so during stable stretches
        # more than half the residuals are exactly zero and the MAD-based
        # robust sigma degenerates to zero -- any nonzero residual would
        # then trivially "exceed" it. Fall back to the plain std, which
        # stays well-defined as long as there's some real scatter (i.e.
        # the dips this is meant to catch).
        sigma = float(np.std(residual))
    return np.abs(residual) > threshold_sigma * sigma


def despike_lightcurve_with_dc(lc: LightCurve, dc_report_path: str, telescope: str = "m2", pixel: int = 0,
                                baseline_window_s: float = 60.0, threshold_sigma: float = 4.0,
                                pad_s: float = 1.0, fill: str = "remove") -> Tuple[LightCurve, float]:
    """
    Use DC excursions (flag_dc_excursions) as a bad-data mask for the fast
    light curve: fast samples within +-pad_s of a flagged DC sample are
    either dropped outright (fill="remove", the default) or linearly
    interpolated over (fill="interpolate").

    Removal is the safer default for *characterising* the noise: an
    interpolated stretch is smooth by construction, so it contributes ~0
    local variance and mechanically pulls down any downstream sigma/PSD
    estimate on the corrected curve -- inflating the apparent improvement
    rather than reflecting a real one. Interpolation fabricates a smooth
    fake signal exactly where a real event could be, too, so a detection
    pipeline should exclude/down-weight these times rather than interpolate
    over them regardless. Removal does leave the time grid non-uniform
    (small gaps where samples were dropped); at the few-percent flagged
    fractions seen so far that's a minor approximation for PSD/binning
    (which both still work off the surviving samples' index order and
    median dt) but would matter more at higher flagged fractions.

    Parameters
    ----------
    lc : LightCurve
        Fast light curve.
    dc_report_path, telescope, pixel : see load_dc_series.
    baseline_window_s, threshold_sigma : see flag_dc_excursions.
    pad_s : float
        Extra time (seconds) padded on each side of a flagged DC sample
        when mapping it onto the fast light curve's time grid.
    fill : {"remove", "interpolate"}
        How to handle flagged fast samples.

    Returns
    -------
    lc_despiked : LightCurve
    bad_fraction : float
        Fraction of fast samples that were flagged.
    """
    if fill not in ("remove", "interpolate"):
        raise ValueError(f"fill must be 'remove' or 'interpolate', got {fill!r}")

    dc_time_s, dc_value = load_dc_series(dc_report_path, telescope=telescope, pixel=pixel)

    # Restrict to this light curve's own observing window before computing the
    # excursion baseline/sigma: the full-night DC series spans long idle/slewing
    # stretches with a completely different (often flat, near-zero-scatter)
    # baseline, which drags the robust sigma down -- in the worst case to exactly
    # zero, which would flag almost every in-window sample as an "excursion".
    dc_dt_s = float(np.median(np.diff(dc_time_s)))
    in_window = (dc_time_s >= lc.time.min() - baseline_window_s) & (dc_time_s <= lc.time.max() + baseline_window_s)
    dc_time_s, dc_value = dc_time_s[in_window], dc_value[in_window]

    bad_dc = flag_dc_excursions(dc_time_s, dc_value, baseline_window_s, threshold_sigma)

    fast_bad = np.zeros(len(lc.time), dtype=bool)
    for t in dc_time_s[bad_dc]:
        fast_bad |= (lc.time >= t - pad_s) & (lc.time <= t + pad_s)

    bad_fraction = float(np.mean(fast_bad))

    meta = dict(lc.meta)
    meta.update({
        "dc_report_path": dc_report_path,
        "dc_telescope": telescope,
        "dc_pixel": pixel,
        "correction": "dc_despike",
        "dc_despike_fill": fill,
        "bad_fraction": bad_fraction,
    })

    if fill == "remove":
        good = ~fast_bad
        return LightCurve(lc.time[good], lc.signal[good], meta=meta), bad_fraction

    signal = lc.signal.copy()
    if np.any(fast_bad) and not np.all(fast_bad):
        good = ~fast_bad
        signal[fast_bad] = np.interp(lc.time[fast_bad], lc.time[good], lc.signal[good])
    return LightCurve(lc.time, signal, meta=meta), bad_fraction
