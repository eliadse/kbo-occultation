# examples/dc_combine_2025_12_16.py
"""
Apply the fast/slow-DC frequency-domain combination (dc_combine.py) to the
2025-12-16 MAGIC-II run and check whether it actually removes the
correlated sub-few-Hz noise excess found in inspect_real_noise_2025_12_16.py.

PSD of the fractional-flux fluctuations, raw vs. DC-combined -- if the
correlated component really is coming from whatever the slow DC already
sees, the combined PSD should show a reduced sub-few-Hz excess.
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

print("\nSaved: noise_psd_dc_combined_2025_12_16.png")
