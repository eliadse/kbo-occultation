# examples/dc_noise_removal_comparison_2025_12_16.py
"""
Compare three ways of using the slow DC report to clean up the fast
2025-12-16 MAGIC-II light curves (see kbo_occultation/dc_combine.py):

  1. fft_swap   -- combine_lightcurve_with_dc: swap in the DC's own
                   frequency content below a cutoff. Already shown (in
                   dc_combine_2025_12_16.py) to make things *worse*, because
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
from kbo_occultation.dc_combine import (
    combine_lightcurve_with_dc,
    detrend_lightcurve_with_dc,
    despike_lightcurve_with_dc,
)

CHANNEL = "A"  # MAGIC-II
DATA_DIR = f"{PACKAGE_DATA}/observations"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"
FREQUENCY_CUT_HZ = 1.0
SMOOTH_WINDOW_S = 10.0
DESPIKE_BASELINE_WINDOW_S = 60.0
DESPIKE_THRESHOLD_SIGMA = 4.0

MAG_B = {
    "HD-4335": 5.95, "HD-8837": 6.55, "HD-22353": 7.04,
    "HD-17316": 7.53, "HD-20920": 7.96,
}


def segmented_welch(time, signal, fs, nperseg_target=2**14, min_run_samples=1024):
    """
    Welch PSD that respects gaps in `time` (e.g. from despike's removed
    samples): splits into contiguous runs (breaks where the local sample
    spacing exceeds 1.5x the nominal dt), Welch's each run separately, and
    averages the results onto a common frequency grid, weighted by run
    length. Falls back to an ordinary single-run Welch when there are no
    gaps.
    """
    dt = 1.0 / fs
    gaps = np.diff(time)
    break_idx = np.where(gaps > 1.5 * dt)[0] + 1
    bounds = np.concatenate([[0], break_idx, [len(time)]])

    freq_common = None
    Pxx_accum = None
    weight_total = 0
    for start, stop in zip(bounds[:-1], bounds[1:]):
        run = signal[start:stop]
        if len(run) < min_run_samples:
            continue
        nperseg = min(nperseg_target, len(run))
        f, Pxx = welch(run, fs=fs, window="hann", nperseg=nperseg, noverlap=nperseg // 2, scaling="density")
        if freq_common is None:
            freq_common = np.linspace(0, fs / 2, nperseg_target // 2 + 1)
        Pxx_interp = np.interp(freq_common, f, Pxx)
        weight = len(run)
        Pxx_accum = Pxx_interp * weight if Pxx_accum is None else Pxx_accum + Pxx_interp * weight
        weight_total += weight

    return freq_common, Pxx_accum / weight_total

filenames = sorted(glob.glob(f"{DATA_DIR}/Spectrum*.bin"))
raw_lcs = {}
for fn in filenames:
    star = re.search(r"15min_(.*)_10002", fn).group(1)
    if star == "dark":
        continue
    raw_lcs[star] = LightCurve.from_stat_binary_all(fn)[CHANNEL]

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
fig1.savefig("noise_psd_comparison_2025_12_16.png", dpi=200)

# ─── Summary table: robust sigma at a few bin sizes ─────────────────────
bin_sizes_report = [1, 100, 5000]
print(f"{'star':10s} {'method':10s}", *[f"n={n:>6d}" for n in bin_sizes_report])
for star in order:
    for label, lc in variants[star].items():
        norm = lc.signal / np.mean(lc.signal)
        row = []
        for n in bin_sizes_report:
            if n == 1:
                row.append(robust_sigma(norm))
            else:
                _, s_b, _ = bin_average(lc.time, norm, n)
                row.append(robust_sigma(s_b))
        print(f"{star:10s} {label:10s}", *[f"{v:10.5f}" for v in row])
    print()

# ─── Noise vs. bin size, all variants overlaid, one panel per star ──────
bin_sizes = np.unique(np.round(np.logspace(0, 4, 30)).astype(int))
fig2, axes2 = plt.subplots(len(order), 1, figsize=(8, 2.6 * len(order)), sharex=True)
for ax, star in zip(np.atleast_1d(axes2), order):
    for label, lc in variants[star].items():
        dt = lc.dt
        norm = lc.signal / np.mean(lc.signal)
        measured, valid_n = [], []
        for n in bin_sizes:
            if n >= len(norm) // 4:
                break
            _, s_b, _ = bin_average(lc.time, norm, n)
            measured.append(robust_sigma(s_b))
            valid_n.append(n)
        valid_n = np.array(valid_n)
        ax.plot(valid_n * dt, measured, lw=1.2, color=colors[label], label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"Scatter $\delta I/I$")
    ax.set_title(star)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
axes2[-1].set_xlabel("Bin duration (s)") if len(order) > 1 else ax.set_xlabel("Bin duration (s)")
fig2.tight_layout()
fig2.savefig("noise_vs_binsize_comparison_2025_12_16.png", dpi=200)

print("Saved: noise_psd_comparison_2025_12_16.png, noise_vs_binsize_comparison_2025_12_16.png")
