"""
Measured signal-to-noise ratio of real fast-photometry observations.

Port of the SNR approaches explored in the old kbo_occultation_rates notebook
(tests/check_fast_photometry_data.ipynb), so that every observation night can
be measured the same way and compared to Instrument.noise_budget:

- signal = mean of the binned light curve, or the last (brightest) peak of
  its histogram (``histogram_peak``), which ignores the wind/tracking dips
  that drag the mean down;
- NSB (+ electronics) level subtracted from it, from a dark run
  (``nsb_from_dark``), from the DC background pixels (``nsb_from_dc_background``)
  or from the first peak of the DC histogram (``nsb_from_dc_histogram``);
- noise = mean per-bin standard error of the binned light curve
  (``measured_snr``, the old "methods 1-3 / p1-p3"), or the robust local
  sigma after removing the slow baseline (``detrended_snr``), which is the
  quantity the theoretical per-sample noise model should be compared to;
- optional good-time windows (``good_time_lightcurve``) and an airmass
  extinction correction (``extinction_factor``).

Superseded pieces of the old notebook that are *not* ported here, because the
package already covers them: ``remove_valleys*`` -> dc_combine.despike_lightcurve_with_dc,
``outlier_removal*`` -> filtering.remove_outliers_lightcurve (applied by
photometry.preprocess_observation), ``combined_dc_reports_with_std_dev_2`` /
``combine_fast_dc_*`` -> dc_combine.combine_lightcurve_with_dc,
``average_samples`` -> detectability.bin_average.
"""

import warnings

import numpy as np
from scipy.signal import welch

from .detectability import bin_average
from .dc_combine import load_dc_series
from .matched_filter import compute_baseline_and_sigma, robust_sigma
from .photometry import LightCurve

# Number of mirror panels per telescope, for the throughput lost to faulty
# AMC panels (throughput_from_faulty_panels). MAGIC: both telescopes use
# magic_spysii's mirror-magic2-nopsf_experimental_clean.def (see
# magic_spysii/reflector.py), 247 facets of 0.985 m side = 239.6 m^2,
# consistent with instruments.mirror_area_per_telescope['MAGIC'] = 240.
N_PANELS = {"m1": 247, "m2": 247, "lst1": 198}


def histogram_peak(data, which="last", n_bins=200, min_count=20):
    """
    Centre of the first local maximum of the histogram of ``data``, scanning
    from the high end (``which="last"``, the star level, below which the
    wind/tracking dips pile up) or from the low end (``which="first"``, e.g.
    the NSB level of a DC series that also contains on-star periods).

    Bins with fewer than ``min_count`` entries are ignored, so isolated
    outliers at either end don't count as a peak. Returns NaN (with a
    warning) when no peak is found.
    """
    counts, edges = np.histogram(np.asarray(data), bins=n_bins)
    centres = 0.5 * (edges[1:] + edges[:-1])
    counts = np.where(counts < min_count, 0, counts)

    order = range(len(counts) - 1, -1, -1) if which == "last" else range(len(counts))
    best_count, best_idx = -1, None
    for i in order:
        if counts[i] >= best_count:
            best_count, best_idx = counts[i], i
        else:
            return float(centres[best_idx])

    warnings.warn("histogram_peak: no peak found; lower min_count")
    return float("nan")


def good_time_lightcurve(lc: LightCurve, window_rel_s):
    """
    Restrict ``lc`` to ``[start, stop]`` seconds after its first sample (the
    observation log's ``good_time``). ``None`` returns ``lc`` unchanged.
    """
    if window_rel_s is None:
        return lc
    t_rel = lc.time - lc.time[0]
    mask = (t_rel >= window_rel_s[0]) & (t_rel <= window_rel_s[1])
    meta = dict(lc.meta, good_time=list(window_rel_s))
    return LightCurve(lc.time[mask], lc.signal[mask], meta=meta)


def nsb_from_dark(lc_dark: LightCurve, n_bin: int):
    """NSB level and its per-bin standard error from a dark run, binned like the stars."""
    _, signal, sigma = bin_average(lc_dark.time, lc_dark.signal, n_bin)
    return float(np.mean(signal)), float(np.mean(sigma))


def nsb_from_dc_background(dc_report_path: str, telescope: str, t_start: float, t_stop: float):
    """
    NSB level (in DC units) from the DC report's background pixels (columns
    1 and 2 for MAGIC; LST-1 only has column 1) during [t_start, t_stop]: the
    pixel with the lower mean is taken, since another may sit on a field star.
    Returns (mean, std).

    Only comparable with signals in DC units, i.e. the output of
    dc_combine.combine_lightcurve_with_dc, not with the raw fast variance.
    """
    levels = []
    for pixel in (1, 2):
        try:
            t, v = load_dc_series(dc_report_path, telescope=telescope, pixel=pixel)
        except IndexError:
            break
        levels.append(v[(t >= t_start) & (t <= t_stop)])
    best = min(levels, key=np.mean)
    return float(np.mean(best)), float(np.std(best))


def nsb_from_dc_histogram(dc_report_path: str, telescope: str, n_bins=100, min_count=20):
    """
    NSB level (in DC units) as the first peak of the histogram of the
    on-target pixel's DC over the whole report -- used for LST-1, whose
    background pixel reads abnormally high. Returns (level, 0.0).
    If the camera was closed between observations, this will return ~0.
    In that case use other method, like the background pixels or a dark patch.
    """
    _, v = load_dc_series(dc_report_path, telescope=telescope, pixel=0)
    return histogram_peak(v, which="first", n_bins=n_bins, min_count=min_count), 0.0


def dc_to_fast_scale(lc: LightCurve, dc_report_path: str, telescope: str) -> float:
    """
    Factor converting DC units into the fast signal's (variance) units for
    this observation: mean(fast) / mean(on-target DC over the same window).
    Lets the DC-based NSB estimates be used on the raw fast signal for nights
    without a dark run. Assumes both are proportional to the photon rate; an
    electronic pedestal in the variance that the DC does not see biases it.
    """
    t, v = load_dc_series(dc_report_path, telescope=telescope, pixel=0)
    sel = (t >= lc.time[0]) & (t <= lc.time[-1])
    return float(np.mean(lc.signal) / np.mean(v[sel]))


def measured_snr(lc: LightCurve, n_bin: int, nsb: float, nsb_err: float = 0.0,
                 signal: str = "mean", peak_kwargs=None) -> dict:
    """
    SNR at the binned cadence (n_bin samples per point):

        SNR = (S - NSB) / sqrt(<sigma_bin>^2 + nsb_err^2)

    with sigma_bin the standard error of each bin (bin_average) and S the mean
    of the binned curve (``signal="mean"``) or the last histogram peak of the
    unbinned one (``signal="peak"``). Returns signal, noise and snr.

    sigma_bin only sees the scatter *within* each bin, so noise correlated on
    the bin timescale or longer is invisible to it (compare detrended_snr);
    and averaging the per-bin std's (ddof=0) underestimates sigma for small
    bins -- by ~8% for n_bin=10 even for pure white noise -- so this SNR runs
    high. Kept for continuity with the old notebook's numbers.
    """
    _, binned, sigma = bin_average(lc.time, lc.signal, n_bin)
    if signal == "mean":
        level = float(np.mean(binned))
    elif signal == "peak":
        level = histogram_peak(lc.signal, which="last", **(peak_kwargs or {}))
    else:
        raise ValueError(f"signal must be 'mean' or 'peak', got {signal!r}")
    star = level - nsb
    noise = float(np.sqrt(np.mean(sigma) ** 2 + nsb_err ** 2))
    return {"signal": star, "noise": noise, "snr": star / noise}


def detrended_snr(lc: LightCurve, n_bin: int, nsb: float, detrend_window_s: float = 10.0) -> dict:
    """
    SNR at the binned cadence from the robust (MAD) bin-to-bin scatter after
    dividing out a rolling-median baseline (the baseline of
    matched_filter.compute_baseline_and_sigma, on the binned curve). Unlike
    measured_snr, drifts slower than the window (transparency, airmass,
    tracking) don't count as noise, while everything faster -- including
    correlated noise on the bin timescale -- does. This is the estimator to
    compare with the per-sample Instrument.noise_budget.

    A single robust sigma over the whole residual is used instead of
    compute_baseline_and_sigma's local sigma: on its decimated path that one
    is a within-block MAD over window/64 samples, blind to noise slower than
    that block, so it would change meaning with the bin size.

    The default 10 s window is ~20x an occultation template, the scale the
    search pipeline detrends on; keep it well below the length of the light
    curve, or it stops removing anything. The baseline-relative sigma is
    converted into a fraction of the star signal (baseline - NSB). Returns
    signal, noise and snr.
    """
    t, binned, _ = bin_average(lc.time, lc.signal, n_bin)
    baseline, _ = compute_baseline_and_sigma(
        t, binned, template_duration_s=detrend_window_s / 20.0,
        detrend_window_s=detrend_window_s)
    sigma_frac = robust_sigma(binned / baseline - 1.0)
    level = float(np.median(baseline))
    star = level - nsb
    noise = sigma_frac * level
    return {"signal": star, "noise": noise, "snr": star / noise}


def airmass_from_zd(zd_deg):
    return 1.0 / np.cos(np.radians(zd_deg))


def extinction_factor(zd_deg, k_mag):
    """
    Factor that brings a signal observed at zenith distance ``zd_deg`` back
    to the zenith, for an extinction coefficient ``k_mag`` (mag / airmass):
    10^(0.4 k (X - 1)).
    """
    return 10 ** (0.4 * k_mag * (airmass_from_zd(zd_deg) - 1.0))


def throughput_from_faulty_panels(n_unfocused, telescope: str) -> float:
    """
    Fraction of mirror area left given the number of panels not correctly
    focused (an observation-log ``faulty_mirrors[telescope]`` value; all
    panels have the same area); 1.0 if unknown (None).
    """
    if n_unfocused is None:
        return 1.0
    return 1.0 - n_unfocused / N_PANELS[telescope]


def band_noise_rms(lc: LightCurve, bands: dict, nperseg: int = 4096) -> dict:
    """
    RMS of the fractional flux (signal / mean - 1) inside frequency bands,
    integrated from its Welch PSD. ``bands`` maps a label to (f_low, f_high) Hz.
    """
    frac = lc.signal / np.mean(lc.signal) - 1.0
    f, pxx = welch(frac, fs=1.0 / lc.dt, window="hann", nperseg=nperseg,
                   noverlap=nperseg // 2, scaling="density")
    trapz = getattr(np, "trapezoid", None) or np.trapz
    out = {}
    for label, (f1, f2) in bands.items():
        sel = (f >= f1) & (f <= f2)
        out[label] = float(np.sqrt(trapz(pxx[sel], f[sel])))
    return out


def noise_vs_bin_size(lc: LightCurve, bin_sizes, nsb: float = 0.0, detrend_window_s: float = 10.0):
    """
    Fractional noise (relative to the star signal) of the binned light curve
    for each bin size, two ways:

    - within: mean per-bin standard error (what measured_snr uses);
    - between: robust bin-to-bin scatter after detrending (1 / detrended_snr).

    For white noise both fall as 1/sqrt(n); where "between" flattens out,
    correlated noise dominates and extra averaging stops helping. Bin sizes
    leaving fewer than 50 bins, or below 2 samples, are skipped.
    Returns (bin_sizes_used, within, between).
    """
    used, within, between = [], [], []
    for n in bin_sizes:
        n = int(n)
        if n < 2 or len(lc.signal) // n < 50:
            continue
        m = measured_snr(lc, n, nsb)
        d = detrended_snr(lc, n, nsb, detrend_window_s=detrend_window_s)
        used.append(n)
        within.append(m["noise"] / m["signal"])
        between.append(1.0 / d["snr"])
    return np.array(used), np.array(within), np.array(between)


def frame_noise(lc: LightCurve) -> float:
    """
    Fractional sample-to-sample scatter of the raw (unbinned) signal,
    robust_sigma(diff(signal)) / sqrt(2) / median(signal). Differencing
    adjacent frames removes everything slower than a couple of frames, so this
    isolates the per-frame noise: photon statistics plus the statistical error
    of the variance estimate itself.
    """
    return robust_sigma(np.diff(lc.signal)) / np.sqrt(2.0) / float(np.median(lc.signal))


def variance_estimator_floor(n_digitizations: int, correlated_samples: float = 1.0) -> float:
    """
    Relative standard error of a variance estimated from ``n_digitizations``
    Gaussian samples, sqrt(2 k / N), with k the number of consecutive samples
    that are correlated (the analog bandwidth spreads a pulse over several
    2 ns samples). It is a fractional floor that does not shrink for brighter
    stars, and the photon-counting noise model does not include it.
    """
    return float(np.sqrt(2.0 * correlated_samples / n_digitizations))
