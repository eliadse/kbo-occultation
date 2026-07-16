# examples/check_highpass_signal_safety.py
"""
Does the 3 Hz high-pass cutoff (filtering.py, see highpass_filter_2025_12_16.py)
clip the KBO occultation signal itself, rather than just the correlated
instrumental/atmospheric noise it was designed to remove?

Two independent checks, for a range of KBO radii spanning what the other
example scripts in this repo already treat as representative
(20-300 m, at 40 AU, 25 km/s shadow velocity -- see plot_diffraction_psd.py,
magic_lst_filter_comparison.py, injection_recovery_test.py):

1. Spectral: what fraction of the noiseless template's own power (in
   temporal frequency, via detectability.spatial_power_spectrum +
   the shadow velocity) sits below the cutoff -- if that's negligible,
   a high-pass filter there can't be removing much of the signal either.

2. Direct: resample the template onto the real fast cadence (~1907 Hz,
   config.standard_sampling/standard_sample_duration), apply
   highpass_butterworth exactly as the real pipeline would, and compare
   the filtered event's matched-filter SNR (detectability.matched_filter_snr,
   scored against the *original* unfiltered template so this measures shape
   preservation, not just amplitude) to the unfiltered value. This is the
   check that actually matters: it captures any distortion from the
   filter's transition band, not just a hard power-below-cutoff cut.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kbo_occultation import (
    GridConfig, NumericalConfig, BandpassConfig, StarConfig, KBOConfig,
    compute_lightcurve,
)
from kbo_occultation.config import standard_sampling, standard_sample_duration
from kbo_occultation.detectability import (
    spatial_to_time, spatial_power_spectrum, resample_to_cadence,
)
from kbo_occultation.filtering import highpass_butterworth

STAR_TEMPERATURE_K = 30000
STAR_ANGULAR_RADIUS_MAS = 0.03
DISTANCE_AU = 40.0
SHADOW_VELOCITY_MPS = 25e3
RADII_M = [20, 50, 80, 100, 150, 200, 300]
CUTOFF_HZ = 3.0
SAMPLE_DT_S = standard_sampling * standard_sample_duration * 1e-9  # real fast-photometry cadence, ~524.3us

BAND = BandpassConfig(lam_min_nm=400, lam_max_nm=500, n_lambda=41)
grid = GridConfig(x_max_m=4000, n_x=4000)
numerics = NumericalConfig(n_r_grid=3000)
star = StarConfig(temperature_K=STAR_TEMPERATURE_K, angular_radius_mas=STAR_ANGULAR_RADIUS_MAS)

print(f"{'R (m)':>6s} {'power<{:.0f}Hz'.format(CUTOFF_HZ):>12s} {'SNR retained':>13s} {'depth orig':>11s} {'depth filt':>11s}")

fig, (ax_time, ax_psd) = plt.subplots(1, 2, figsize=(13, 5.5))

for radius_m in RADII_M:
    kbo = KBOConfig(radius_m=radius_m, distance_au=DISTANCE_AU)
    x_m, intensity = compute_lightcurve(kbo, star, BAND, grid, numerics)

    # --- 1. spectral check: fraction of power below the cutoff ---
    freq_per_m, power = spatial_power_spectrum(x_m, intensity)
    freq_hz = freq_per_m * SHADOW_VELOCITY_MPS
    cumulative = np.cumsum(power) / np.sum(power)
    below_cut = float(np.interp(CUTOFF_HZ, freq_hz, cumulative))

    # --- 2. direct check: resample to real cadence, filter, compare ---
    t_s = spatial_to_time(x_m, SHADOW_VELOCITY_MPS)
    t_grid, I_grid = resample_to_cadence(t_s, intensity, SAMPLE_DT_S)

    filtered = highpass_butterworth(I_grid, fs=1.0 / SAMPLE_DT_S, cutoff_hz=CUTOFF_HZ)

    # Score the *filtered* curve's shape against the *original* template's
    # own residual -- this is what actually tells us whether filtering
    # distorted the event, vs. just restating its own amplitude.
    residual_orig = I_grid - 1.0
    residual_filt = filtered - 1.0
    snr_ratio = float(np.sum(residual_filt * residual_orig) / np.sum(residual_orig**2))

    depth_orig = 1.0 - I_grid.min()
    depth_filt = 1.0 - filtered.min()

    print(f"{radius_m:6d} {100*below_cut:11.4f}% {100*snr_ratio:12.2f}% {depth_orig:11.4f} {depth_filt:11.4f}")

    ax_psd.plot(freq_hz, power, lw=1.2, label=f"R={radius_m}m")

    if radius_m in (max(RADII_M), 80):
        i_min = np.argmin(I_grid)
        mask = slice(max(0, i_min - 2000), min(len(t_grid), i_min + 2000))
        ax_time.plot(1e3 * (t_grid[mask] - t_grid[i_min]), I_grid[mask], lw=1.3, label=f"R={radius_m}m orig")
        ax_time.plot(1e3 * (t_grid[mask] - t_grid[i_min]), filtered[mask], "--", lw=1.3, label=f"R={radius_m}m filtered")

ax_psd.axvline(CUTOFF_HZ, color="k", ls="--", lw=1.0, alpha=0.7, label="cutoff")
ax_psd.set_xscale("log")
ax_psd.set_yscale("log")
ax_psd.set_xlim(1e-2, 3 * SHADOW_VELOCITY_MPS / (2 * grid.x_max_m / grid.n_x))
ax_psd.set_xlabel("Temporal frequency (Hz)")
ax_psd.set_ylabel("Normalised PSD")
ax_psd.set_title("Template PSD vs. cutoff")
ax_psd.legend(fontsize=7)
ax_psd.grid(alpha=0.3, which="both")

ax_time.set_xlabel("Time from dip minimum (ms)")
ax_time.set_ylabel("Relative intensity")
ax_time.set_title("Original vs. high-pass-filtered template (worst + fiducial case)")
ax_time.legend(fontsize=8)
ax_time.grid(alpha=0.3)

fig.tight_layout()
fig.savefig("highpass_signal_safety_check.png", dpi=200)
print("\nSaved: highpass_signal_safety_check.png")
