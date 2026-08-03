# examples/highpass_filter_2025_12_16.py
"""
High-pass filtering (kbo_occultation/filtering.py) of the 2025-12-16
MAGIC-II light curves, independent of the DC-based corrections in
dc_combine.py: it only needs a uniform time grid and a cutoff frequency.

Cutoff: inspect_real_noise_2025_12_16.py found the correlated noise excess
rolling over into the white floor around 3-5 Hz; CUTOFF_HZ=3.0 sits at the
conservative (lower) end of that so as to not eat into the occultation
signal's own band above it.

Three variants per star:
  1. raw
  2. highpass       -- highpass_lightcurve(raw, CUTOFF_HZ)
  3. highpass+despike -- despike_lightcurve_with_dc applied *after* the
     high-pass filter (not before: despike's default fill="remove" leaves
     gaps, and a high-pass filter assumes uniform sampling like any
     FFT/IIR method would -- see filtering.py's module docstring).

This demonstrates the two corrections compose (in that order) rather than
picking one over the other.
"""

import glob
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import welch

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import LightCurve
from kbo_occultation.filtering import highpass_lightcurve
from kbo_occultation.dc_combine import despike_lightcurve_with_dc
from kbo_occultation.detectability import bin_average
from kbo_occultation.matched_filter import compute_baseline_and_sigma

CHANNEL = "A"  # MAGIC-II
DATA_DIR = f"{PACKAGE_DATA}/observations"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"
CUTOFF_HZ = 3.5
DISPLAY_STRIDE = 20  # thin the ~1.7M-point series for legible time-domain plots
BLOCK_SAMPLES = 200  # samples per point in the binned mean+-error overlay

MAG_B = {
    "HD-4335": 5.95, "HD-8837": 6.55, "HD-22353": 7.04,
    "HD-17316": 7.53, "HD-20920": 7.96,
}


def segmented_welch(time, signal, fs, nperseg_target=2**14, min_run_samples=1024):
    """See dc_noise_removal_comparison_2025_12_16.py -- gap-aware Welch,
    needed once despike (fill="remove") is in the chain."""
    dt = 1.0 / fs
    gaps = np.diff(time)
    break_idx = np.where(gaps > 1.5 * dt)[0] + 1
    bounds = np.concatenate([[0], break_idx, [len(time)]])

    freq_common, Pxx_accum, weight_total = None, None, 0
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

filenames = sorted(glob.glob(f"{DATA_DIR}/Spectrum*20251216T*.npz"))
raw_lcs = {}
for fn in filenames:
    star = re.search(r"15min_(.*)_10002", fn).group(1)
    if star == "dark":
        continue
    raw_lcs[star] = LightCurve.from_preprocessed_all(fn)[CHANNEL]

order = sorted(raw_lcs, key=lambda s: MAG_B.get(s, 99))

variants = {}
for star in order:
    lc = raw_lcs[star]
    hp = highpass_lightcurve(lc, CUTOFF_HZ)
    hp_despike, bad_frac = despike_lightcurve_with_dc(hp, DC_PKL)
    variants[star] = {"raw": lc, "highpass": hp, "highpass+despike": hp_despike}

# ─── PSD grid ──────────────────────────────────────────────────────────
colors = {"raw": "C0", "highpass": "C1", "highpass+despike": "C3"}
fig1, axes1 = plt.subplots(len(order), 1, figsize=(9, 2.6 * len(order)), sharex=True)
for ax, star in zip(np.atleast_1d(axes1), order):
    for label, lc in variants[star].items():
        fs = 1.0 / lc.dt
        norm = lc.signal / np.mean(raw_lcs[star].signal) - 1.0
        f, Pxx = segmented_welch(lc.time, norm, fs)
        ax.plot(f[1:], Pxx[1:], lw=1.1, color=colors[label], label=label)
    ax.axvline(CUTOFF_HZ, color="k", ls="--", lw=1.0, alpha=0.6, label="cutoff")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(r"PSD of $\delta I/I$")
    ax.set_title(f"{star} (B={MAG_B[star]})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
axes1[-1].set_xlabel("Frequency (Hz)") if len(order) > 1 else ax.set_xlabel("Frequency (Hz)")
fig1.tight_layout()
fig1.savefig("noise_psd_highpass_2025_12_16.png", dpi=200)

# ─── Time-domain grid ──────────────────────────────────────────────────
# Left: raw counts (each star at its own absolute level, drifting over the
# run). Right: the same series after the high-pass
sl = slice(None, None, DISPLAY_STRIDE)
fig2, axes2 = plt.subplots(len(order), 2, figsize=(14, 2.4 * len(order)), sharex=True)
axes2 = np.atleast_2d(axes2)
for row, star in zip(axes2, order):
    ax_raw, ax_hp = row
    lc = variants[star]["raw"]
    hp = variants[star]["highpass"]
    hpd = variants[star]["highpass+despike"]

    t0 = lc.time.min()

    ax_raw.plot(lc.time[sl] - t0, lc.signal[sl], ".", ms=0.8, color=colors["raw"])
    tb, sb, eb = bin_average(lc.time, lc.signal, BLOCK_SAMPLES)
    #base, sig = compute_baseline_and_sigma(lc.time, lc.signal, 0.01, 0.01)
    # TODO use this method for the sigma that you want to compare, it makes no sense to 
    # create another function when you want to compare the bin_average error with the 
    # error from the compute_baseline_and_sigma. Change this in the following curves as well
    ax_raw.errorbar(tb - t0, sb, yerr=robust_bin_sem(lc.signal, BLOCK_SAMPLES),
                    fmt="none", ecolor="C2", elinewidth=1.2, capsize=2, zorder=4,
                    label=r"robust (MAD)/$\sqrt{N}$")
    ax_raw.errorbar(tb - t0, sb, yerr=eb, fmt=".", ms=2.5, color="k",
                    ecolor="k", elinewidth=0.6, capsize=1, zorder=3,
                    label=f"{BLOCK_SAMPLES}-sample mean, std/$\\sqrt{{N}}$")
    ax_raw.set_ylabel(f"{star}\n(B={MAG_B[star]})\nraw counts")
    ax_raw.grid(alpha=0.3)
    ax_raw.legend(fontsize=7, loc="upper right", markerscale=4)

    ax_hp.plot(hp.time[sl] - t0, (hp.signal)[sl],
               ".", ms=0.8, color=colors["highpass"], label="highpass")
    ax_hp.plot(hpd.time[sl] - t0, (hpd.signal)[sl],
               ".", ms=0.8, color=colors["highpass+despike"], label="highpass+despike")
    
    tb, sb, eb = bin_average(hp.time, hp.signal, BLOCK_SAMPLES)

    ax_hp.errorbar(tb - t0, sb, yerr=robust_bin_sem(hp.signal, BLOCK_SAMPLES),
                   fmt="none", ecolor="C2", elinewidth=1.2, capsize=2, zorder=4,
                   label=r"robust (MAD)/$\sqrt{N}$")
    ax_hp.errorbar(tb - t0, sb, yerr=eb, fmt=".", ms=2.5, color="k",
                   ecolor="k", elinewidth=0.6, capsize=1, zorder=3,
                   label=f"{BLOCK_SAMPLES}-sample mean, std/$\\sqrt{{N}}$")
    ax_hp.axhline(0, color="k", lw=0.6)
    ax_hp.set_ylabel(r"$I$")
    # robust y-limits so a few despike spikes don't flatten the rest
    resid = hp.signal
    ax_hp.set_ylim(*np.percentile(resid, [0.1, 99.9]))
    ax_hp.grid(alpha=0.3)
    ax_hp.legend(fontsize=7, loc="upper right", markerscale=6)

axes2[0, 0].set_title("Raw light curve")
axes2[0, 1].set_title(f"After high-pass ({CUTOFF_HZ} Hz) -- slow drift removed")
axes2[-1, 0].set_xlabel("Time (s)")
axes2[-1, 1].set_xlabel("Time (s)")
fig2.tight_layout()
fig2.savefig("highpass_timedomain_2025_12_16.png", dpi=140)

plt.show()

print(f"Saved: noise_psd_highpass_2025_12_16.png (cutoff={CUTOFF_HZ} Hz)")
