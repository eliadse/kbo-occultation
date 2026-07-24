# examples/baseline_fit_test_2025_12_16.py
"""
Visual sanity check of matched_filter.compute_baseline_and_sigma on the real
2025_12_16 observations.

compute_baseline_and_sigma is the workhorse of the matched-filter search: it
fits a rolling-median baseline to a light curve and reports the robust
per-sample *fractional* noise sigma. It is also the most expensive step in the
search and uses a decimated-median approximation to stay fast. This script
checks, by eye, that the baseline it produces actually tracks each real light
curve and that the reported sigma is a sensible noise level -- by overlaying,
per star (channel A = MAGIC-II, the only real channel that night):

    * the raw variance light curve,
    * the fitted baseline, and
    * a +/-1 sigma band.

Note the band is drawn as baseline * (1 +/- sigma), NOT baseline +/- sigma:
sigma is the robust std of the *fractional* residual (signal/baseline - 1), so
in signal units the 1-sigma spread scales with the baseline.
"""

import glob
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import LightCurve
from kbo_occultation.matched_filter import compute_baseline_and_sigma

DATA_DIR = f"{PACKAGE_DATA}/observations"
CHANNEL = "A"  # MAGIC-II, the only real channel on 2025_12_16

# Rolling-median baseline window. Fixed (rather than derived from a specific KBO
# template) so it is consistent across stars and comfortably longer than any
# real sub-second event -- wide enough not to distort an event, per
# compute_baseline_and_sigma's docstring.
DETREND_WINDOW_S = 2.0
TEMPLATE_DURATION_S = 0.1  # nominal; only floors the median window when detrend_window_s is given

DISPLAY_STRIDE = 10  # decimate for plotting only; baseline/sigma use full resolution

# Brightest first, dark last (magnitudes from inspect_real_noise_2025_12_16.py).
MAG_B = {
    "HD-4335": 5.95, "HD-8837": 6.55, "HD-22353": 7.04,
    "HD-17316": 7.53, "HD-20920": 7.96,
}

# ─── Load the 2025_12_16 light curves (the 20251216T timestamp excludes the
#     2026 run files that share this directory) ──────────────────────────────
filenames = sorted(glob.glob(f"{DATA_DIR}/Spectrum*20251216T*.npz"))
lightcurves = {}
for fn in filenames:
    star = re.search(r"15min_(.*)_10002", fn).group(1)
    lightcurves[star] = LightCurve.from_preprocessed_all(fn)[CHANNEL]

order = sorted(lightcurves, key=lambda s: MAG_B.get(s, 99))  # dark (not in MAG_B) sorts last

# ─── Fit the baseline for each star ──────────────────────────────────────────
fits = {}
print(f"{'star':10s} {'n':>9s} {'sigma':>10s} {'noise%':>8s}")
for star in order:
    lc = lightcurves[star]
    baseline, sigma = compute_baseline_and_sigma(
        lc.time, lc.signal, TEMPLATE_DURATION_S, detrend_window_s=DETREND_WINDOW_S
    )
    fits[star] = (baseline, sigma)
    print(f"{star:10s} {len(lc.time):9d} {np.median(sigma):10.5f}")

# ─── Plot: one panel per star ────────────────────────────────────────────────
fig, axes = plt.subplots(len(order), 1, figsize=(12, 3 * len(order)), sharex=True)
sl = slice(None, None, DISPLAY_STRIDE)

for ax, star in zip(np.atleast_1d(axes), order):
    lc = lightcurves[star]
    baseline, sigma = fits[star]
    # Relative time: the stars were recorded sequentially through the night, so
    # their absolute timestamps don't overlap -- rezero each to its own start so
    # every panel shows its full ~426 s span across the full width.
    t = lc.time - lc.time.min()

    ax.plot(t[sl], lc.signal[sl], ".", ms=1, alpha=0.8, color="C0", label="raw")
    ax.fill_between(
        t[sl], (baseline * (1 - sigma))[sl], (baseline * (1 + sigma))[sl],
        color="C1", alpha=0.2, label="±1σ",
    )
    ax.plot(t[sl], baseline[sl], color="C1", lw=1.2, label="baseline")

    ax.set_title(f"{star}  (σ={np.median(sigma):.3g})")
    ax.set_ylabel("Signal")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

np.atleast_1d(axes)[-1].set_xlabel("Time (s)")
fig.tight_layout()
fig.savefig("baseline_fit_test_2025_12_16.png", dpi=300)
plt.show()
print("\nSaved baseline_fit_test_2025_12_16.png")
