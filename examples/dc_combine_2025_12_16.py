# examples/dc_combine_2025_12_16.py
"""
Apply the fast/slow-DC frequency-domain combination (dc_combine.py) to the
2025-12-16 MAGIC-II run and check whether it actually removes the
correlated sub-few-Hz noise excess found in inspect_real_noise_2025_12_16.py.

Two panels per star:
  1. PSD of the fractional-flux fluctuations, raw vs. DC-combined.
  2. Noise-vs-bin-size curve, raw vs. DC-combined, against the white-noise
     1/sqrt(n) prediction -- if the correlated component really is coming
     from whatever the slow DC already sees, the combined curve should sit
     much closer to that prediction over the range where averaging
     previously did nothing.
"""

import glob
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import LightCurve
from kbo_occultation.detectability import bin_average
from kbo_occultation.matched_filter import robust_sigma
from kbo_occultation.dc_combine import combine_lightcurve_with_dc

CHANNEL = "A"  # MAGIC-II
DATA_DIR = f"{PACKAGE_DATA}/observations"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"
FREQUENCY_CUT_HZ = 1.0

MAG_B = {
    "HD-4335": 5.95, "HD-8837": 6.55, "HD-22353": 7.04,
    "HD-17316": 7.53, "HD-20920": 7.96,
}

filenames = sorted(glob.glob(f"{DATA_DIR}/Spectrum*.npz"))
raw_lcs = {}
for fn in filenames:
    star = re.search(r"15min_(.*)_10002", fn).group(1)
    lc = LightCurve.from_preprocessed_all(fn)[CHANNEL]
    if star == "dark":
        continue  # no star flux -> DC report has nothing meaningful to combine
    raw_lcs[star] = lc

order = sorted(raw_lcs, key=lambda s: MAG_B.get(s, 99))

combined_lcs = {}
for star in order:
    combined_lcs[star] = combine_lightcurve_with_dc(
        raw_lcs[star], DC_PKL, frequency_cut_hz=FREQUENCY_CUT_HZ, telescope="m2", pixel=0
    )

# ─── PSD: raw vs. combined ──────────────────────────────────────────────
fig1, axes1 = plt.subplots(len(order), 1, figsize=(9, 2.6 * len(order)), sharex=True)
for ax, star in zip(np.atleast_1d(axes1), order):
    for label, lc, style in [("raw", raw_lcs[star], "-"), ("DC-combined", combined_lcs[star], "-")]:
        fs = 1.0 / lc.dt
        norm = lc.signal / np.mean(lc.signal) - 1.0
        nperseg = min(2**16, len(norm))
        f, Pxx = welch(norm, fs=fs, window="hann", nperseg=nperseg, noverlap=nperseg // 2, scaling="density")
        ax.plot(f[1:], Pxx[1:], style, lw=1.1, label=label)
    ax.axvline(FREQUENCY_CUT_HZ, color="k", ls="--", lw=1.0, alpha=0.6, label="cut")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"PSD of $\delta I/I$")
    ax.set_title(f"{star} (B={MAG_B[star]})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
axes1[-1].set_xlabel("Frequency (Hz)") if len(order) > 1 else ax.set_xlabel("Frequency (Hz)")
fig1.tight_layout()
fig1.savefig("noise_psd_dc_combined_2025_12_16.png", dpi=200)

# ─── Noise vs. bin size: raw vs. combined ──────────────────────────────
fig2, ax2 = plt.subplots(figsize=(8, 6))
bin_sizes = np.unique(np.round(np.logspace(0, 4, 30)).astype(int))

print(f"{'star':10s} {'raw n=1':>10s} {'raw n=5000':>12s} {'comb n=1':>10s} {'comb n=5000':>12s}")
for star in order:
    for label, lc, ls in [("raw", raw_lcs[star], "--"), ("DC-combined", combined_lcs[star], "-")]:
        dt = lc.dt
        norm = lc.signal / np.mean(lc.signal)
        sigma1 = robust_sigma(norm)

        measured, valid_n = [], []
        for n in bin_sizes:
            if n >= len(norm) // 4:
                break
            _, s_b, _ = bin_average(lc.time, norm, n)
            measured.append(robust_sigma(s_b))
            valid_n.append(n)
        valid_n = np.array(valid_n)

        color = f"C{order.index(star)}"
        ax2.plot(valid_n * dt, measured, ls, lw=1.3, color=color,
                  label=f"{star} ({label})" if label == "raw" else None,
                  alpha=1.0 if label == "DC-combined" else 0.4)

        if label == "raw":
            row = [sigma1]
        else:
            row.append(sigma1)

    # print n=1 / n=5000 for both raw and combined
    raw_sigma1 = robust_sigma(raw_lcs[star].signal / np.mean(raw_lcs[star].signal))
    comb_sigma1 = robust_sigma(combined_lcs[star].signal / np.mean(combined_lcs[star].signal))
    _, raw_5000, _ = bin_average(raw_lcs[star].time, raw_lcs[star].signal / np.mean(raw_lcs[star].signal), 5000)
    _, comb_5000, _ = bin_average(combined_lcs[star].time, combined_lcs[star].signal / np.mean(combined_lcs[star].signal), 5000)
    print(f"{star:10s} {raw_sigma1:10.5f} {robust_sigma(raw_5000):12.5f} {comb_sigma1:10.5f} {robust_sigma(comb_5000):12.5f}")

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("Bin duration (s)")
ax2.set_ylabel(r"Scatter of binned $\delta I/I$ (robust $\sigma$)")
ax2.set_title("Noise vs. averaging timescale: raw (dashed, faint) vs. DC-combined (solid)")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3, which="both")
fig2.tight_layout()
fig2.savefig("noise_vs_binsize_dc_combined_2025_12_16.png", dpi=200)

print("\nSaved: noise_psd_dc_combined_2025_12_16.png, noise_vs_binsize_dc_combined_2025_12_16.png")
