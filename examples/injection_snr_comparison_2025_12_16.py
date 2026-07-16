# examples/injection_snr_comparison_2025_12_16.py
"""
Does either noise-reduction approach developed for the 2025-12-16 MAGIC-II
run (kbo_occultation/filtering.py's high-pass, or kbo_occultation/dc_combine.py's
DC-based corrections) actually improve recovery of a real occultation signal,
rather than just reducing noise in the abstract?

Injects a simulated KBO template (as in tests/injection_recovery_test.py) into
the raw HD-17316 light curve, applies each correction to the *injected*
series (the order a real pipeline would use -- correct first, search after),
and runs the same blind matched-filter search
(matched_filter.sliding_matched_filter_snr) on every variant. Comparing the
recovered SNR at the known truth time isolates what each correction actually
does to a real signal, rather than just to the noise floor in the abstract.

Variants (see dc_noise_removal_comparison_2025_12_16.py and
highpass_filter_2025_12_16.py for the PSD-level picture behind each):
  raw               -- no correction.
  highpass          -- highpass_lightcurve, cutoff below the occultation
                       band (see highpass_filter_2025_12_16.py).
  dc_fft_swap       -- combine_lightcurve_with_dc; already shown to make the
                       noise floor *worse* below ~1Hz (the DC's own sample
                       scatter exceeds the fast signal's) -- kept in as a
                       check that a worse noise floor really does cost SNR
                       here too, not just look worse in a PSD.
  dc_detrend        -- detrend_lightcurve_with_dc.
  dc_despike        -- despike_lightcurve_with_dc (fill="remove"); leaves
                       gaps, so this variant's false-alarm count away from
                       the injection is only approximate (see dc_combine.py's
                       module docstring on segmented PSDs) -- the recovery
                       peak itself is unaffected as long as no despiked
                       sample falls inside the injection window (checked
                       below).
  highpass_despike  -- despike applied *after* highpass, per filtering.py's
                       ordering note.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import LightCurve
from kbo_occultation.instruments import Instrument
from kbo_occultation import (
    GridConfig, NumericalConfig, BandpassConfig, StarConfig, KBOConfig,
    compute_lightcurve, inject_occultation, sliding_matched_filter_snr, find_candidates,
)
from kbo_occultation.detectability import spatial_to_time
from kbo_occultation.filtering import highpass_lightcurve
from kbo_occultation.dc_combine import (
    combine_lightcurve_with_dc, detrend_lightcurve_with_dc, despike_lightcurve_with_dc,
)

DATA_DIR = f"{PACKAGE_DATA}/observations"
DATA_FILE = f"{DATA_DIR}/Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_HD-17316_10002_20002_20251216T223319.bin"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"
CHANNEL = "A"

STAR_TEMPERATURE_K = 30000
STAR_ANGULAR_RADIUS_MAS = 0.03
DISTANCE_AU = 40.0
SHADOW_VELOCITY_MPS = 25e3
RADIUS_M = 80
SITE = "La Palma"
SNR_THRESHOLD = 5.0
# Reduced-chi2 cutoff for the shape veto (matched_filter.shape_veto_chi2).
MAX_CHI2_REDUCED = 2.0

HIGHPASS_CUTOFF_HZ = 3.0  # matches highpass_filter_2025_12_16.py
DC_FREQUENCY_CUT_HZ = 1.0  # matches dc_combine_2025_12_16.py
DETREND_SMOOTH_WINDOW_S = 10.0
DESPIKE_BASELINE_WINDOW_S = 60.0
DESPIKE_THRESHOLD_SIGMA = 4.0

grid = GridConfig(x_max_m=4000, n_x=800)
numerics = NumericalConfig(n_r_grid=3000)

# ─── 1. Simulate the noiseless template ────────────────────────────────
inst = Instrument(telescope_type="MAGIC", filter_type=None, site=SITE)
l_min, l_max, n_lambda = inst.response_band()
band = BandpassConfig(l_min, l_max, n_lambda)

kbo = KBOConfig(radius_m=RADIUS_M, distance_au=DISTANCE_AU)
star = StarConfig(temperature_K=STAR_TEMPERATURE_K, angular_radius_mas=STAR_ANGULAR_RADIUS_MAS)

x_m, intensity = compute_lightcurve(kbo, star, band, grid, numerics, response=inst.total_transmission)
template_time = spatial_to_time(x_m, SHADOW_VELOCITY_MPS)
template_duration_s = float(template_time.max() - template_time.min())

print(f"Template: R={RADIUS_M}m depth={1 - intensity.min():.4f} duration={template_duration_s:.4f}s")

# ─── 2. Load real data, inject once at a fixed, known t0 ───────────────
lc = LightCurve.from_stat_binary_all(DATA_FILE)[CHANNEL]
t0_s = lc.time[0] + 0.5 * (lc.time[-1] - lc.time[0])
inj = inject_occultation(lc.time, lc.signal, x_m, intensity, SHADOW_VELOCITY_MPS, t0_s=t0_s)
injected_lc = LightCurve(inj.time, inj.injected_signal, meta=dict(lc.meta))

print(f"Injected at t0={inj.t0_s:.3f}s, window={inj.window_s}")

# ─── 3. Apply each correction to the injected series ───────────────────
hp = highpass_lightcurve(injected_lc, HIGHPASS_CUTOFF_HZ)
fft_swap = combine_lightcurve_with_dc(injected_lc, DC_PKL, frequency_cut_hz=DC_FREQUENCY_CUT_HZ)
detrended = detrend_lightcurve_with_dc(injected_lc, DC_PKL, smooth_window_s=DETREND_SMOOTH_WINDOW_S)
despiked, bad_frac = despike_lightcurve_with_dc(
    injected_lc, DC_PKL, baseline_window_s=DESPIKE_BASELINE_WINDOW_S, threshold_sigma=DESPIKE_THRESHOLD_SIGMA
)
hp_despiked, hp_bad_frac = despike_lightcurve_with_dc(
    hp, DC_PKL, baseline_window_s=DESPIKE_BASELINE_WINDOW_S, threshold_sigma=DESPIKE_THRESHOLD_SIGMA
)

variants = {
    "raw": injected_lc,
    "highpass": hp,
    "dc_fft_swap": fft_swap,
    "dc_detrend": detrended,
    "dc_despike": despiked,
    "highpass_despike": hp_despiked,
}

for label, bf in (("dc_despike", bad_frac), ("highpass_despike", hp_bad_frac)):
    in_window_before = np.sum((lc.time >= inj.window_s[0]) & (lc.time <= inj.window_s[1]))
    in_window_after = np.sum((variants[label].time >= inj.window_s[0]) & (variants[label].time <= inj.window_s[1]))
    print(f"{label}: flagged {100*bf:.2f}% of all samples, "
          f"removed {in_window_before - in_window_after} of {in_window_before} inside the injection window")

# ─── 4. Same blind matched-filter search for every variant ─────────────
print(f"\n{'method':18s} {'sigma':>10s} {'recovered SNR':>14s} {'offset (ms)':>12s} "
      f"{'candidates':>11s} {'false alarms':>13s}")

results = {}
for label, v_lc in variants.items():
    mf = sliding_matched_filter_snr(v_lc.time, v_lc.signal, template_time, intensity)
    candidates = find_candidates(mf, snr_threshold=SNR_THRESHOLD)
    survivors = [c for c in candidates if c.chi2_reduced <= MAX_CHI2_REDUCED]
    results[label] = (mf, candidates, survivors)

    if survivors:
        nearest = min(survivors, key=lambda c: abs(c.time_s - inj.t0_s))
        offset_ms = (nearest.time_s - inj.t0_s) * 1e3
        false_alarms = len(survivors) - 1
        print(f"{label:18s} {mf.sigma:10.5f} {nearest.snr:14.2f} {offset_ms:12.3f} "
              f"{len(candidates):11d} {false_alarms:13d}")
    else:
        print(f"{label:18s} {mf.sigma:10.5f} {'none':>14s} {'--':>12s} {len(candidates):11d} {'--':>13s}")

# ─── Plots ──────────────────────────────────────────────────────────────
pad_s = 5 * template_duration_s
fig, axes = plt.subplots(len(variants), 1, figsize=(9, 2.2 * len(variants)), sharex=True)
for ax, (label, v_lc) in zip(axes, variants.items()):
    mask = (v_lc.time > inj.t0_s - pad_s) & (v_lc.time < inj.t0_s + pad_s)
    ax.plot(v_lc.time[mask], v_lc.signal[mask], ".", ms=2, alpha=0.4, color="C0")
    ax.axvline(inj.t0_s, color="k", ls="--", lw=1)
    _, _, survivors = results[label]
    if survivors:
        nearest = min(survivors, key=lambda c: abs(c.time_s - inj.t0_s))
        ax.set_title(f"{label} (recovered SNR={nearest.snr:.2f})")
    else:
        ax.set_title(f"{label} (not recovered)")
    ax.set_ylabel("Signal")
axes[-1].set_xlabel("Time (s)")
fig.suptitle(f"Injected R={RADIUS_M}m event, zoomed on truth window, per noise-reduction variant")
fig.tight_layout()
fig.savefig("injection_snr_comparison_2025_12_16.png", dpi=200)

fig2, ax2 = plt.subplots(figsize=(7, 4))
labels = list(variants.keys())
snrs = []
for label in labels:
    _, _, survivors = results[label]
    if survivors:
        nearest = min(survivors, key=lambda c: abs(c.time_s - inj.t0_s))
        snrs.append(nearest.snr)
    else:
        snrs.append(0.0)
ax2.bar(labels, snrs, color="C0")
ax2.axhline(SNR_THRESHOLD, color="gray", ls="--", lw=1, label=f"threshold={SNR_THRESHOLD}")
ax2.set_ylabel("Recovered matched-filter SNR")
ax2.set_title(f"Injected R={RADIUS_M}m event: recovered SNR by noise-reduction method")
plt.setp(ax2.get_xticklabels(), rotation=30, ha="right")
ax2.legend(fontsize=8)
fig2.tight_layout()
fig2.savefig("injection_snr_bar_2025_12_16.png", dpi=200)

print("\nSaved: injection_snr_comparison_2025_12_16.png, injection_snr_bar_2025_12_16.png")
