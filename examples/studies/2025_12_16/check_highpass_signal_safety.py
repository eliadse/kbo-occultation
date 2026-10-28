# examples/studies/2025_12_16/check_highpass_signal_safety.py
"""
Does the high-pass (config.HIGHPASS_CUTOFF_HZ) clip the occultation signal itself?

The template is embedded in a long flat record on the real cadence before
filtering (a cutoff period comparable to the array being filtered gives
meaningless results, which is what the earlier version of this check did),
the filter is applied exactly as in the pipeline, and the event is scored
against the *original* template. Run over a scan of cutoffs and radii:

- power fraction: share of the template's power below the cutoff (spectral view);
- amplitude retained: <filtered, original> / <original, original>;
- SNR retained: matched-filter SNR of the filtered event against the
  *reshaped* template (what the pipeline's search actually uses) divided by
  the unfiltered SNR, for white noise. 1 means the filter costs nothing.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kbo_occultation import GridConfig, NumericalConfig, BandpassConfig, StarConfig, KBOConfig, compute_lightcurve
from kbo_occultation.config import HIGHPASS_CUTOFF_HZ, standard_sampling, standard_sample_duration
from kbo_occultation.analysis.injection import inject_occultation
from kbo_occultation.observations.filtering import highpass_butterworth

OUT_DIR = os.path.join("outputs", "studies/2025_12_16")
os.makedirs(OUT_DIR, exist_ok=True)
STAR = StarConfig(temperature_K=10000, angular_radius_mas=0.03)
DISTANCE_AU = 39.0
V_MPS = 25e3
RADII_M = [20, 50, 100, 200, 300, 1000]
CUTOFFS_HZ = [0.3, 0.5, 1.0, 2.0, HIGHPASS_CUTOFF_HZ, 5.0]
DT = standard_sampling * standard_sample_duration * 1e-9
RECORD_S = 240.0                     # >> 1 / lowest cutoff
BAND = BandpassConfig(320.0, 450.0, 21)
GRID = GridConfig(x_max_m=10000, n_x=4000)
NUM = NumericalConfig(n_r_grid=3000, n_star_side=48)

time = np.arange(int(RECORD_S / DT)) * DT
t0 = time[len(time) // 2]
flat = np.ones_like(time)
print(f"{'R (m)':>6s} {'cut Hz':>7s} {'power<cut':>10s} {'amp kept':>9s} {'SNR kept':>9s}")
snr_kept = np.zeros((len(RADII_M), len(CUTOFFS_HZ)))
for i, radius in enumerate(RADII_M):
    x_m, intensity = compute_lightcurve(KBOConfig(radius, DISTANCE_AU), STAR, BAND, GRID, NUM)
    inj = inject_occultation(time, flat, x_m, intensity, V_MPS, t0_s=t0)
    orig = inj.injected_signal - 1.0
    spec = np.abs(np.fft.rfft(orig)) ** 2          # on the embedded record: resolution 1/RECORD_S
    freq = np.fft.rfftfreq(len(orig), DT)
    cum = np.cumsum(spec) / np.sum(spec)
    for j, fc in enumerate(CUTOFFS_HZ):
        filt = highpass_butterworth(inj.injected_signal, 1.0 / DT, fc) - np.mean(highpass_butterworth(flat, 1.0 / DT, fc))
        shape = filt - np.mean(filt[:1000])                    # reshaped template (flat record has none of the noise)
        below = float(np.interp(fc, freq, cum))
        amp = float(np.sum(shape * orig) / np.sum(orig ** 2))
        # white-noise matched filter on the reshaped template: SNR ~ |shape|; unfiltered ~ |orig|
        snr = float(np.sqrt(np.sum(shape ** 2) / np.sum(orig ** 2)))
        snr_kept[i, j] = snr
        print(f"{radius:6d} {fc:7.2f} {100*below:9.2f}% {100*amp:8.1f}% {100*snr:8.1f}%")

fig, ax = plt.subplots(figsize=(7, 4.5))
for i, radius in enumerate(RADII_M):
    ax.plot(CUTOFFS_HZ, snr_kept[i], "o-", label=f"R={radius} m")
ax.axvline(HIGHPASS_CUTOFF_HZ, color="k", ls="--", lw=1, label="pipeline cutoff")
ax.set_xscale("log"); ax.set_xlabel("high-pass cutoff [Hz]"); ax.set_ylabel("white-noise SNR retained")
ax.set_ylim(0, 1.05); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/highpass_signal_safety_check.png", dpi=200)
print("Saved: highpass_signal_safety_check.png")
