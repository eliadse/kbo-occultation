# examples/studies/2025_12_16/dc_noise_removal_comparison.py
"""
Compare three ways of using the slow DC report to clean up the fast
2025-12-16 MAGIC-II light curves (see kbo_occultation/dc_combine.py):

  1. fft_swap   -- combine_lightcurve_with_dc: swap in the DC's own
                   frequency content below a cutoff. Already shown (in
                   examples/studies/2025_12_16/dc_combine.py) to make things *worse*, because
                   the DC's raw per-sample scatter is larger than the fast
                   signal's own low-frequency scatter.
  2. detrend    -- detrend_lightcurve_with_dc: divide by a heavily-smoothed
                   (median-filtered) DC trend instead of the raw values, so
                   only the genuine multi-second drift is imported, not the
                   DC's own sample noise.
  3. despike    -- despike_lightcurve_with_dc: use excursions in the DC as
                   a bad-data mask. Flagged samples are *removed* (not
                   interpolated over): an interpolated stretch is smooth by
                   construction and would mechanically pull down any
                   downstream sigma/PSD estimate, inflating the apparent
                   improvement rather than reflecting a real one.

Removal leaves gaps in the time grid, though -- naively concatenating the
surviving samples and handing them to a single Welch call splices unrelated
segments together, and the resulting phase discontinuities spray spurious
power across a broad range of higher frequencies (checked directly: with
despike's ~11 gaps spread through ~900k samples, most 65536-sample Welch
segments span at least one splice). segmented_welch below instead runs
Welch separately on each contiguous run between gaps and averages the
result, weighted by run length -- for the raw/fft_swap/detrend curves
(no gaps) this is identical to a plain Welch call.

For each star: PSD of all four variants (raw + the three corrections), and
a summary table of robust-sigma noise at a few bin sizes.
"""

import os
import glob
import re

import numpy as np
from kbo_occultation.analysis.spectra import segmented_welch as _segmented_welch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.observations.photometry import LightCurve
from kbo_occultation.observations.dc_combine import (
    combine_lightcurve_with_dc,
    detrend_lightcurve_with_dc,
    despike_lightcurve_with_dc,
)

OUT_DIR = os.path.join("outputs", "studies/2025_12_16")
os.makedirs(OUT_DIR, exist_ok=True)
CHANNEL = "A"  # MAGIC-II
DATA_DIR = f"{PACKAGE_DATA}/observations"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"
from kbo_occultation.config import DC_CUT_HZ, DC_DETREND_WINDOW_S
FREQUENCY_CUT_HZ = DC_CUT_HZ["m2"]
SMOOTH_WINDOW_S = DC_DETREND_WINDOW_S
DESPIKE_BASELINE_WINDOW_S = 60.0
DESPIKE_THRESHOLD_SIGMA = 4.0

MAG_B = {
    "HD-4335": 5.95, "HD-8837": 6.55, "HD-22353": 7.04,
    "HD-17316": 7.53, "HD-20920": 7.96,
}


def segmented_welch(time, signal, fs, nperseg_target=2**14, min_run_samples=1024):
    # moved to kbo_occultation.analysis.spectra (default segment there is 2**16); this study keeps 2**14.
    return _segmented_welch(time, signal, fs, nperseg_target, min_run_samples)

filenames = sorted(glob.glob(f"{DATA_DIR}/Spectrum*.npz"))
raw_lcs = {}
for fn in filenames:
    star = re.search(r"15min_(.*)_10002", fn).group(1)
    if star == "dark":
        continue
    raw_lcs[star] = LightCurve.from_preprocessed_all(fn)[CHANNEL]

order = sorted(raw_lcs, key=lambda s: MAG_B.get(s, 99))

variants = {}
bad_fractions = {}
for star in order:
    lc = raw_lcs[star]
    fft_swap = combine_lightcurve_with_dc(lc, DC_PKL, frequency_cut_hz=FREQUENCY_CUT_HZ)
    detrended = detrend_lightcurve_with_dc(lc, DC_PKL, smooth_window_s=SMOOTH_WINDOW_S)
    despiked, bad_frac = despike_lightcurve_with_dc(
        lc, DC_PKL, baseline_window_s=DESPIKE_BASELINE_WINDOW_S, threshold_sigma=DESPIKE_THRESHOLD_SIGMA
    )
    variants[star] = {"raw": lc, "fft_swap": fft_swap, "detrend": detrended, "despike": despiked}
    bad_fractions[star] = bad_frac

# ─── PSD grid ──────────────────────────────────────────────────────────
colors = {"raw": "C0", "fft_swap": "C1", "detrend": "C2", "despike": "C3"}
fig1, axes1 = plt.subplots(len(order), 1, figsize=(9, 2.6 * len(order)), sharex=True)
for ax, star in zip(np.atleast_1d(axes1), order):
    for label, lc in variants[star].items():
        fs = 1.0 / lc.dt
        norm = lc.signal / np.mean(lc.signal) - 1.0
        f, Pxx = segmented_welch(lc.time, norm, fs)
        ax.plot(f[1:], Pxx[1:], lw=1.1, color=colors[label], label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"PSD of $\delta I/I$")
    ax.set_title(f"{star} (B={MAG_B[star]}, despike flagged {100*bad_fractions[star]:.2f}% of samples)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
axes1[-1].set_xlabel("Frequency (Hz)") if len(order) > 1 else ax.set_xlabel("Frequency (Hz)")
fig1.tight_layout()
fig1.savefig(f"{OUT_DIR}/noise_psd_comparison_2025_12_16.png", dpi=200)

print("Saved: noise_psd_comparison_2025_12_16.png")
