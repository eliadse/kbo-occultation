# kbo_occultation/detectability.py
"""
Turn a simulated diffraction light curve into detectability metrics.

These functions are deliberately noise-model-agnostic: they take a per-sample noise sigma 
(in the same normalised-intensity units as the light curve, i.e. fractional flux) as input 
rather than computing it internally. 

That's so they can plug into whatever noise characterization you settle on 
-- the three-component model (photon noise, read noise, Young's-formula scintillation, 
systematic floor combined in quadrature) discussed for the injection pipeline, or just a 
flat sigma for quick tests -- without this module needing to know about it.
"""

from typing import Optional, Tuple

import numpy as np


def spatial_to_time(x_m: np.ndarray, shadow_velocity_mps: float) -> np.ndarray:
    """
    Convert the spatial diffraction-pattern axis to time, given the KBO
    shadow's speed across the observer's location. t=0 at x=0 (closest
    approach to the shadow center).

    Typical TNO shadow velocities are ~20-25 km/s.
    """
    return np.asarray(x_m) / shadow_velocity_mps


def spatial_power_spectrum(x_m: np.ndarray, intensity: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Power spectrum of a diffraction pattern I(x), evaluated on its native
    uniform spatial grid.

    Returns
    -------
    freq_cyc_per_m : ndarray
        Non-negative spatial frequencies (cycles/m).
    power : ndarray
        Power spectral density, normalised so its peak is 1.
    """
    x_m = np.asarray(x_m)
    dx = x_m[1] - x_m[0]

    signal = intensity - np.mean(intensity)
    spectrum = np.abs(np.fft.rfft(signal)) ** 2
    freq_cyc_per_m = np.fft.rfftfreq(len(x_m), d=dx)

    return freq_cyc_per_m, spectrum / spectrum.max()


def nyquist_frequency(x_m: np.ndarray, intensity: np.ndarray, shadow_velocity_mps: float,
                       power_threshold: float = 0.99) -> Tuple[float, float]:
    """
    Estimate the Nyquist sampling frequency needed to resolve a simulated
    diffraction light curve, from the pattern's own spatial frequency
    content -- rather than assuming the Fresnel scale alone sets it.

    Must be called on the *spatial* pattern I(x_m), on its native uniform
    grid, before any resampling to a fixed time cadence (resample_to_cadence
    can only remove high-frequency content, never add it back, so measuring
    bandwidth after resampling just tells you the cadence you already
    chose, not what the physics needs).

    The pattern's spatial bandwidth (highest spatial frequency f_x with
    power above `power_threshold` of the peak) is converted to a temporal
    frequency via f_t = f_x * shadow_velocity_mps, since x = v*t for a
    shadow crossing the observer at constant velocity.

    Parameters
    ----------
    x_m : array_like
        Uniform spatial grid (m).
    intensity : array_like
        Noiseless diffraction pattern I(x_m).
    shadow_velocity_mps : float
        KBO shadow velocity across the observer's location (m/s).
    power_threshold : float
        Fraction of peak spectral power above which a spatial frequency
        counts as "present" in the signal. Default 1%.

    Returns
    -------
    f_nyquist_hz : float
        Highest temporal frequency present in the pattern. Sample at
        >= 2 * f_nyquist_hz to avoid aliasing (in practice, use a few x
        margin -- power_threshold is a soft cutoff, not a hard band limit).
    dt_max_s : float
        Corresponding maximum sampling interval, 1 / (2 * f_nyquist_hz).
    """
    freq_cyc_per_m, power = spatial_power_spectrum(x_m, intensity)
    cumulative = np.cumsum(power) / np.sum(power)
    above = freq_cyc_per_m[cumulative > power_threshold]
    f_max_spatial = float(above[0]) if above.size else 0.0

    f_nyquist_hz = f_max_spatial * shadow_velocity_mps
    dt_max_s = 1.0 / (2.0 * f_nyquist_hz) if f_nyquist_hz > 0 else np.inf

    return f_nyquist_hz, dt_max_s


def resample_to_cadence(t_s: np.ndarray, intensity: np.ndarray, dt_s: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Resample a light curve onto a fixed sampling cadence dt_s (e.g. an
    instrument's real integration time) by linear interpolation.
    """
    t_s = np.asarray(t_s)
    t_new = np.arange(t_s.min(), t_s.max(), dt_s)
    I_new = np.interp(t_new, t_s, intensity)
    return t_new, I_new


def bin_average(t_s: np.ndarray, signal: np.ndarray, n_samples: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Bin a *noisy* light curve by averaging non-overlapping chunks of
    n_samples raw points (any leftover points that don't fill a final
    chunk are dropped) -- unlike resample_to_cadence's interpolation,
    this actually reduces per-sample noise (by sqrt(n_samples), to the
    extent the noise is white on the timescale being binned) and cuts
    the number of independent samples a downstream search has to try,
    which is the point when raw cadence is much finer than the features
    you're actually looking for.

    Use nyquist_frequency on the noiseless template to find how coarse
    you can safely go without smearing out the feature itself (its
    dt_max_s), rather than picking n_samples by guesswork.

    Parameters
    ----------
    t_s, signal : array_like
        Uniformly-sampled light curve.
    n_samples : int
        Number of consecutive raw samples averaged into each output
        point.

    Returns
    -------
    t_binned, signal_binned, sigma_binned : ndarray
        sigma_binned is the standard error of the mean within each bin
        (``std(chunk) / sqrt(n_samples)``) -- an empirical, per-bin
        noise estimate.
    """
    t_s = np.asarray(t_s)
    signal = np.asarray(signal)

    n_trim = len(signal) - (len(signal) % n_samples)
    t_chunks = t_s[:n_trim].reshape(-1, n_samples)
    signal_chunks = signal[:n_trim].reshape(-1, n_samples)

    t_binned = t_chunks.mean(axis=1)
    signal_binned = signal_chunks.mean(axis=1)
    sigma_binned = signal_chunks.std(axis=1) / np.sqrt(n_samples)

    return t_binned, signal_binned, sigma_binned


def peak_snr(intensity: np.ndarray, sigma: float) -> float:
    """
    Simplest detectability statistic: the deepest single-sample dip,
    divided by the per-sample noise sigma.

    Fast and conservative, but ignores signal spread over neighbouring
    samples -- for small/fast KBOs where the dip only spans a couple of
    samples this is close to optimal; for wider dips, matched_filter_snr
    makes better use of the data.
    """
    depth = 1.0 - np.min(intensity)
    return float(depth / sigma)


def matched_filter_snr(intensity: np.ndarray, sigma: float, template: Optional[np.ndarray] = None):
    """
    Integrated (matched-filter) SNR of the whole dip, combining
    information from every sample instead of just the deepest one:

        SNR = sqrt( sum_i (1 - template_i)^2 ) / sigma

    assuming uniform, uncorrelated per-sample noise sigma. 

    If `template` is omitted, `intensity` is used as its own template (appropriate
    when testing against the noiseless simulated curve itself, e.g. to map out 
    detectability across a parameter grid before touching real data).
    """
    ref = np.asarray(intensity if template is None else template)
    residual = ref - 1.0

    return float(np.sqrt(np.sum(residual**2)) / sigma)


def event_duration(x_m: np.ndarray, intensity: np.ndarray, threshold: float = 0.5) -> float:
    """
    Full width, in the same units as x_m, over which the intensity drops below `threshold` 
    fraction of the maximum depth (default: full width at half depth). Returns 0.0 if the 
    dip never reaches that threshold.
    """
    x_m = np.asarray(x_m)
    intensity = np.asarray(intensity)
    depth = 1.0 - np.min(intensity)
    if depth <= 0:
        return 0.0

    level = 1.0 - threshold * depth
    below = intensity < level
    if not np.any(below):
        return 0.0

    idx = np.where(below)[0]

    return float(x_m[idx[-1]] - x_m[idx[0]])


def is_detectable(intensity: np.ndarray, sigma: float, snr_threshold: float = 5.0, method: str = "matched_filter") -> bool:
    """
    Convenience wrapper returning a boolean detection flag.

    method : {"matched_filter", "peak"}
    """
    if method == "matched_filter":
        snr = matched_filter_snr(intensity, sigma)
    elif method == "peak":
        snr = peak_snr(intensity, sigma)
    else:
        raise ValueError(f"Unknown method: {method!r}")

    return snr >= snr_threshold


def sigma_from_instrument(instrument, magnitude: float, time_binning, ntels: int = 1) -> float:
    """
    Convert an Instrument's photon-counting SNR (instruments.Instrument.signal_to_noise_ratio) 
    into a per-sample fractional-flux noise sigma, for use as the `sigma` argument to
    peak_snr / matched_filter_snr / is_detectable.

    sigma = 1 / SNR_photometric, where
    SNR_photometric = star_photons / sqrt(star_photons + NSB_photons)
    over `time_binning`.

    Caveat: this currently captures photon shot noise (star + NSB) only.
    It does NOT yet include read noise, scintillation (Young's formula), or a systematic noise floor 
    -- treat resulting SNR/detectability numbers as optimistic upper bounds, not final answers, 
    until the full noise model is wired in.

    Parameters
    ----------
    instrument : instruments.Instrument
        Configured with the telescope_type / filter_type / site you want to evaluate. 
        instrument.signal_to_noise_ratio must be able to run
        (i.e. its required reference data files must be present).
    magnitude : float
        Apparent magnitude of the target star, in the instrument's current filter_type 
        (or unfiltered response, if filter_type is None).
    time_binning : astropy.units.Quantity
        Sampling cadence, e.g. ``524.288 * u.us``.
    ntels : int
        Number of telescopes combined (as in Instrument.signal_to_noise_ratio).
    """
    snr_phot = instrument.signal_to_noise_ratio(magnitude, time_binning, ntels=ntels)
    return float(1.0 / snr_phot)
