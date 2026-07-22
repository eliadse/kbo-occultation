# examples/inspect_real_noise_2025_12_16.py
"""
Noise characterization of the real 2025-12-16 fast-photometry light curves
(MAGIC-II, channel A -- the only telescope with real data that night).

Previously we tried binning/averaging the raw ~1907 Hz signal to cut the
number of matched-filter trials and improve SNR (see
detectability.bin_average, tests/injection_recovery_test.py), implicitly
assuming the per-sample noise is white so the scatter shrinks as
1/sqrt(n_samples). This script checks that assumption directly against
the real data:

1. PSD of each star's (and dark's) fractional-flux fluctuations -- a flat
   PSD is white noise; excess power at low frequency (red/1/f) or narrow
   spectral lines (periodic contamination) both violate the white-noise
   assumption and would explain why averaging under-delivers.
2. A noise-vs-bin-size ("Allan-like") curve: the actual scatter of
   bin-averaged points vs. bin size, against the 1/sqrt(n) white-noise
   prediction. Where the two curves diverge marks the timescale beyond
   which extra averaging stops helping.
3. Coherence between the fast MAGIC-II signal (binned down to ~1 Hz) and
   the slow DC report (pixel 1, the on-target pixel) over the same
   window, to check whether the low-frequency excess is explained by
   whatever is already visible in the slow monitoring (wind/tracking/
   cloud transparency changes).
"""

import glob
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import welch, coherence

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import LightCurve
from kbo_occultation.detectability import bin_average
from kbo_occultation.matched_filter import robust_sigma

CHANNEL = "A"  # MAGIC-II
DATA_DIR = f"{PACKAGE_DATA}/observations"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"

MAG_B = {
    "HD-4335": 5.95, "HD-8837": 6.55, "HD-22353": 7.04,
    "HD-17316": 7.53, "HD-20920": 7.96,
}

# ─── Load the 6 light curves ──────────────────────────────────────────
filenames = sorted(glob.glob(f"{DATA_DIR}/Spectrum*.bin"))
lightcurves = {}
for fn in filenames:
    star = re.search(r"15min_(.*)_10002", fn).group(1)
    lightcurves[star] = LightCurve.from_stat_binary_all(fn)[CHANNEL]

order = sorted(lightcurves, key=lambda s: MAG_B.get(s, 99))  # brightest first, dark last

print(f"{'star':10s} {'n':>8s} {'dt_ms':>8s} {'mean':>10s} {'std':>10s} {'noise%':>8s}")
for star in order:
    lc = lightcurves[star]
    dt = lc.dt
    mean, std = np.mean(lc.signal), np.std(lc.signal)
    print(f"{star:10s} {len(lc.signal):8d} {dt*1e3:8.4f} {mean:10.3f} {std:10.3f} {100*std/mean:8.3f}")

# ─── 1. PSD of fractional-flux fluctuations ───────────────────────────
fig1, ax1 = plt.subplots(figsize=(9, 6))
for star in order:
    lc = lightcurves[star]
    fs = 1.0 / lc.dt
    norm = lc.signal / np.mean(lc.signal) - 1.0
    nperseg = min(2**16, len(norm))
    f, Pxx = welch(norm, fs=fs, window="hann", nperseg=nperseg, noverlap=nperseg // 2, scaling="density")
    ax1.plot(f[1:], Pxx[1:], lw=1.1, label=f"{star} (B={MAG_B[star]})" if star in MAG_B else star)

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel(r"PSD of $\delta I / I$ (Hz$^{-1}$)")
ax1.set_title("2025-12-16, MAGIC-II (channel A): fractional-flux PSD")
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3, which="both")
fig1.tight_layout()
fig1.savefig("noise_psd_2025_12_16.png", dpi=200)

# ─── 2. Noise vs. bin size (Allan-like curve) ──────────────────────────
fig2, ax2 = plt.subplots(figsize=(8, 6))
bin_sizes = np.unique(np.round(np.logspace(0, 4, 30)).astype(int))
bin_sizes = bin_sizes[bin_sizes >= 1]

for star in order:
    lc = lightcurves[star]
    dt = lc.dt
    norm = lc.signal / np.mean(lc.signal)
    sigma1 = robust_sigma(norm)  # single-sample scatter, robust to outliers/glitches

    measured = []
    valid_n = []
    for n in bin_sizes:
        if n >= len(norm) // 4:
            break
        _, s_b, _ = bin_average(lc.time, norm, n)
        measured.append(robust_sigma(s_b))
        valid_n.append(n)

    valid_n = np.array(valid_n)
    line, = ax2.plot(valid_n * dt, measured, marker=".", ms=4, lw=1, label=star)
    ax2.plot(valid_n * dt, sigma1 / np.sqrt(valid_n), ls="--", lw=0.8, color=line.get_color(), alpha=0.6)

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("Bin duration (s)")
ax2.set_ylabel(r"Scatter of binned $\delta I/I$ (robust $\sigma$)")
ax2.set_title("Noise vs. averaging timescale (dashed = white-noise $1/\\sqrt{n}$ prediction)")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3, which="both")
fig2.tight_layout()
fig2.savefig("noise_vs_binsize_2025_12_16.png", dpi=200)

# ─── 3. Coherence between fast signal and slow DC ──────────────────────
try:
    import pickle
    from astropy.time import Time

    with open(DC_PKL, "rb") as fh:
        dc_report = pickle.load(fh)

    dc_times = Time(dc_report.m2["times"], format="isot", scale="utc")
    dc_values = np.array(dc_report.m2["dcs"][:, 0])  # on-target pixel

    fig3, axes3 = plt.subplots(len(order), 1, figsize=(9, 2.2 * len(order)), sharex=False)
    for ax, star in zip(np.atleast_1d(axes3), order):
        lc = lightcurves[star]
        t_fast_unix = lc.time  # already unix seconds (from time_stamp/1e6)
        t0 = Time(t_fast_unix[0], format="unix", scale="utc")
        t1 = Time(t_fast_unix[-1], format="unix", scale="utc")

        mask = (dc_times >= t0) & (dc_times <= t1)
        if mask.sum() < 8:
            ax.set_title(f"{star}: not enough overlapping DC samples")
            continue

        # Bin the fast signal down to the DC's ~1 Hz cadence for a fair comparison
        dt_fast = lc.dt
        n_bin = max(1, int(round(1.0 / dt_fast)))
        t_binned, s_binned, _ = bin_average(lc.time, lc.signal, n_bin)

        dc_sel_times = dc_times[mask].unix
        dc_sel_values = dc_values[mask]
        dc_on_fast = np.interp(t_binned, dc_sel_times, dc_sel_values,
                                left=np.nan, right=np.nan)
        good = ~np.isnan(dc_on_fast)

        fs_binned = 1.0 / (t_binned[1] - t_binned[0])
        nperseg = min(256, good.sum())
        if nperseg < 8:
            ax.set_title(f"{star}: not enough binned samples for coherence")
            continue
        f_coh, Cxy = coherence(s_binned[good], dc_on_fast[good], fs=fs_binned, nperseg=nperseg)

        ax.plot(f_coh, Cxy, lw=1.2)
        ax.set_ylim(-0.03, 1.03)  # so Cxy~1 doesn't sit invisibly on the frame
        ax.set_ylabel("Coherence")
        ax.set_title(f"{star}: fast (binned to {fs_binned:.2f} Hz) vs. slow DC (pixel 1)")
        ax.grid(alpha=0.3)
    axes3[-1].set_xlabel("Frequency (Hz)") if len(order) > 1 else ax.set_xlabel("Frequency (Hz)")
    fig3.tight_layout()
    fig3.savefig("noise_dc_coherence_2025_12_16.png", dpi=200)
except Exception as exc:  # pragma: no cover - diagnostic script
    print(f"Skipping DC coherence check: {exc!r}")

print("\nSaved: noise_psd_2025_12_16.png, noise_vs_binsize_2025_12_16.png, noise_dc_coherence_2025_12_16.png")
