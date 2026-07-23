# examples/injection_montecarlo_2025_12_16.py
"""
Monte Carlo extension of injection_snr_comparison_2025_12_16.py: instead of one 
KBO radius and one random injection time, sweep RADII_M x N_TRIALS_PER_RADIUS 
random t0s (same shared t0 list reused across every radius and method, so 
comparisons stay paired) for each of the 6 noise-reduction variants, and look at 
the resulting SNR/recovery-rate distributions instead of one noisy point estimate.

The single-injection script's result was inconclusive on its own: recovered SNR 
didn't track each method's known noise quality at all (dc_fft_swap, independently 
known to have a *worse* noise floor, scored higher than highpass, which has the 
best) -- a single-realization artifact. The one credible finding was in false-alarm 
counts (a whole-series statistic, not tied to one instant): highpass cut them ~9x, 
the DC-based corrections didn't. This script gets proper statistics on both.

Performance: see kbo_occultation/injection_batch.py's module docstring. In short, 
a full blind search costs ~93s (dominated by a rolling-median baseline), so this 
batch reuses each method's baseline/sigma -- computed once from that method's full, 
*uninjected* light curve -- as an approximation for the noise floor of every 
trial's local scoring window, rather than recomputing it per trial. This is a known 
approximation (a short injected dip could, in principle, shift a median-filtered 
baseline computed directly on top of it), not an exact computation. The validation
spot-check below checks how much this matters in practice, before the full grid 
runs: it compares this fast, cached-baseline local path against a slow, full-series 
recompute (no shortcuts at all) for one trial.

Caveat worth keeping in mind reading the results: the real HD-17316 recording is 
only ~426s long, and each trial needs a pad_s=20s margin, so the ~20 trials per 
radius are not fully independent draws -- they overlap and re-examine the same 
handful of real noise features. Treat the resulting distributions as "how recovery 
varies with when in this run the event lands," not a rigorous i.i.d. ensemble.
"""

import time as time_module

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import LightCurve
from kbo_occultation.instruments import Instrument
from kbo_occultation import (
    GridConfig, NumericalConfig, BandpassConfig, StarConfig, KBOConfig,
    compute_lightcurve, run_injection_monte_carlo,
)
from kbo_occultation.injection_batch import to_dataframe, build_search_template
from kbo_occultation.injection import inject_occultation, random_injection_time
from kbo_occultation.detectability import spatial_to_time
from kbo_occultation.matched_filter import sliding_matched_filter_snr, find_candidates, compute_baseline_and_sigma
from kbo_occultation.filtering import highpass_lightcurve
from kbo_occultation.dc_combine import (
    combine_lightcurve_with_dc, detrend_lightcurve_with_dc, despike_lightcurve_with_dc,
)

DATA_DIR = f"{PACKAGE_DATA}/observations"
DATA_FILE = f"{DATA_DIR}/Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_HD-17316_10002_20002_20251216T223319.npz"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"
CHANNEL = "A"

STAR_TEMPERATURE_K = 30000
STAR_ANGULAR_RADIUS_MAS = 0.03
DISTANCE_AU = 40.0
SHADOW_VELOCITY_MPS = 25e3
SITE = "La Palma"
SNR_THRESHOLD = 5.0
MAX_CHI2_REDUCED = 2.0

HIGHPASS_CUTOFF_HZ = 3.0
DC_FREQUENCY_CUT_HZ = 1.0
DETREND_SMOOTH_WINDOW_S = 10.0
DESPIKE_BASELINE_WINDOW_S = 60.0
DESPIKE_THRESHOLD_SIGMA = 4.0

RADII_M = [50, 100, 150, 200, 300]
N_TRIALS_PER_RADIUS = 20  # change this to run more/fewer trials
PAD_S = 20.0
SEED = 1216

grid = GridConfig(x_max_m=4000, n_x=800)
numerics = NumericalConfig(n_r_grid=3000)

inst = Instrument(telescope_type="MAGIC", filter_type=None, site=SITE)
l_min, l_max, n_lambda = inst.response_band()
band = BandpassConfig(l_min, l_max, n_lambda)
star = StarConfig(temperature_K=STAR_TEMPERATURE_K, angular_radius_mas=STAR_ANGULAR_RADIUS_MAS)

method_specs = {
    "raw": lambda lc: lc,
    "highpass": lambda lc: highpass_lightcurve(lc, HIGHPASS_CUTOFF_HZ),
    "dc_fft_swap": lambda lc: combine_lightcurve_with_dc(lc, DC_PKL, frequency_cut_hz=DC_FREQUENCY_CUT_HZ),
    "dc_detrend": lambda lc: detrend_lightcurve_with_dc(lc, DC_PKL, smooth_window_s=DETREND_SMOOTH_WINDOW_S),
    "dc_despike": lambda lc: despike_lightcurve_with_dc(
        lc, DC_PKL, baseline_window_s=DESPIKE_BASELINE_WINDOW_S, threshold_sigma=DESPIKE_THRESHOLD_SIGMA
    )[0],
    "highpass_despike": lambda lc: despike_lightcurve_with_dc(
        highpass_lightcurve(lc, HIGHPASS_CUTOFF_HZ), DC_PKL,
        baseline_window_s=DESPIKE_BASELINE_WINDOW_S, threshold_sigma=DESPIKE_THRESHOLD_SIGMA
    )[0],
}
METHODS = list(method_specs.keys())

lc = LightCurve.from_preprocessed_all(DATA_FILE)[CHANNEL]

# ─── Validation spot-check: cached-baseline local path vs. slow full-series ground truth ───
print("Validation spot-check (1 trial, method='highpass', radius=RADII_M[0]) ...")
val_radius = RADII_M[0]
kbo = KBOConfig(radius_m=val_radius, distance_au=DISTANCE_AU)
x_m, intensity = compute_lightcurve(kbo, star, band, grid, numerics, response=inst.total_transmission)
template_time = spatial_to_time(x_m, SHADOW_VELOCITY_MPS)
template_duration_s = float(template_time.max() - template_time.min())

val_rng = np.random.default_rng(SEED)
val_t0 = random_injection_time(lc.time, x_m, SHADOW_VELOCITY_MPS, margin_s=PAD_S, rng=val_rng)
inj = inject_occultation(lc.time, lc.signal, x_m, intensity, SHADOW_VELOCITY_MPS, t0_s=val_t0)
injected_lc = LightCurve(inj.time, inj.injected_signal, meta=dict(lc.meta))
corrected = method_specs["highpass"](injected_lc)

# Same per-method "search template" the batch itself now uses (see
# injection_batch.build_search_template) -- matches highpass's own reshaping of the
# event instead of always comparing against the raw, unfiltered template.
template_pad_s = 10.0 * template_duration_s
search_time, search_intensity = build_search_template(
    lc, method_specs["highpass"], x_m, intensity, SHADOW_VELOCITY_MPS, template_pad_s
)

t_start = time_module.perf_counter()
# detrend_window_s must be set from the *original* event duration, not search_time's own
# (deliberately wider) span -- otherwise the default (20x whatever template span it's given)
# blows up into a huge median-filter window.
mf_full = sliding_matched_filter_snr(
    corrected.time, corrected.signal, search_time, search_intensity,
    detrend_window_s=20.0 * template_duration_s,
)
full_time_s = time_module.perf_counter() - t_start
survivors_full = [c for c in find_candidates(mf_full, snr_threshold=SNR_THRESHOLD) if c.chi2_reduced <= MAX_CHI2_REDUCED]

uninjected_corrected = method_specs["highpass"](lc)
cached_baseline, cached_sigma = compute_baseline_and_sigma(
    uninjected_corrected.time, uninjected_corrected.signal, template_duration_s
)
mask = (corrected.time >= val_t0 - PAD_S) & (corrected.time <= val_t0 + PAD_S)
t_start = time_module.perf_counter()
mf_local = sliding_matched_filter_snr(
    corrected.time[mask], corrected.signal[mask], search_time, search_intensity,
    baseline=cached_baseline[mask], sigma=cached_sigma,
)
local_time_s = time_module.perf_counter() - t_start
survivors_local = [c for c in find_candidates(mf_local, snr_threshold=SNR_THRESHOLD) if c.chi2_reduced <= MAX_CHI2_REDUCED]

nearest_full = min(survivors_full, key=lambda c: abs(c.time_s - val_t0)) if survivors_full else None
nearest_local = min(survivors_local, key=lambda c: abs(c.time_s - val_t0)) if survivors_local else None
if nearest_full and nearest_local:
    rel_diff = abs(nearest_full.snr - nearest_local.snr) / nearest_full.snr
    print(f"  full-series path:  SNR={nearest_full.snr:.3f}  ({full_time_s:.1f}s)")
    print(f"  cached-local path: SNR={nearest_local.snr:.3f}  ({local_time_s:.1f}s)")
    print(f"  relative difference: {100 * rel_diff:.2f}%\n")
else:
    print(f"  full-series recovered: {bool(nearest_full)}, cached-local recovered: {bool(nearest_local)}\n")

# ─── Full Monte Carlo grid ────────────────────────────────────────────────
print(f"Running Monte Carlo: {len(RADII_M)} radii x {N_TRIALS_PER_RADIUS} trials x {len(METHODS)} methods ...")
t_start = time_module.perf_counter()
trial_results, noise_profiles = run_injection_monte_carlo(
    lc=lc,
    radii_m=RADII_M,
    n_trials_per_radius=N_TRIALS_PER_RADIUS,
    grid=grid,
    numerics=numerics,
    star=star,
    band=band,
    response=inst.total_transmission,
    shadow_velocity_mps=SHADOW_VELOCITY_MPS,
    distance_au=DISTANCE_AU,
    method_specs=method_specs,
    pad_s=PAD_S,
    snr_threshold=SNR_THRESHOLD,
    max_chi2_reduced=MAX_CHI2_REDUCED,
    seed=SEED,
    checkpoint_path="injection_montecarlo_trials_2025_12_16.csv",
    reshape_template_methods=["highpass", "highpass_despike"],
    per_trial_reshape_methods=["dc_fft_swap"],
)
print(f"Done in {time_module.perf_counter() - t_start:.1f}s")

trial_df = to_dataframe(trial_results)
noise_df = to_dataframe(noise_profiles)
trial_df.to_csv("injection_montecarlo_trials_2025_12_16.csv", index=False)
noise_df.to_csv("injection_montecarlo_noise_profiles_2025_12_16.csv", index=False)

# ─── Plot 1: recovered-SNR distribution, box per (radius, method) ─────────
n_methods = len(METHODS)
width = 0.8 / n_methods
colors = plt.cm.tab10(np.linspace(0, 1, n_methods))
x = np.arange(len(RADII_M))

fig1, ax1 = plt.subplots(figsize=(11, 6))
for i, method in enumerate(METHODS):
    data, positions = [], []
    for j, radius in enumerate(RADII_M):
        vals = trial_df[(trial_df.method == method) & (trial_df.radius_m == radius)]["recovered_snr"].dropna().values
        data.append(vals if len(vals) else [np.nan])
        positions.append(j + (i - (n_methods - 1) / 2) * width)
    bp = ax1.boxplot(data, positions=positions, widths=width * 0.9, patch_artist=True, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(colors[i])
        patch.set_alpha(0.6)
ax1.set_xticks(x)
ax1.set_xticklabels(RADII_M)
ax1.axhline(SNR_THRESHOLD, color="gray", ls="--", lw=1)
ax1.set_xlabel("KBO radius (m)")
ax1.set_ylabel("Recovered matched-filter SNR")
ax1.set_title(f"Recovered SNR distribution ({N_TRIALS_PER_RADIUS} trials/radius, NaN dropped)")
legend_handles = [plt.Rectangle((0, 0), 1, 1, fc=colors[i], alpha=0.6) for i in range(n_methods)]
ax1.legend(legend_handles, METHODS, fontsize=8, loc="upper left")
fig1.tight_layout()
fig1.savefig("injection_montecarlo_snr_boxplot_2025_12_16.png", dpi=200)

# ─── Plot 2: recovery rate (fraction with recovered=True) per (radius, method) ──
fig2, ax2 = plt.subplots(figsize=(9, 5))
for i, method in enumerate(METHODS):
    rates = []
    for radius in RADII_M:
        sub = trial_df[(trial_df.method == method) & (trial_df.radius_m == radius)]
        rates.append(sub["recovered"].mean() if len(sub) else np.nan)
    ax2.bar(x + (i - (n_methods - 1) / 2) * width, rates, width=width * 0.9, color=colors[i], label=method)
ax2.set_xticks(x)
ax2.set_xticklabels(RADII_M)
ax2.set_ylim(0, 1.05)
ax2.set_xlabel("KBO radius (m)")
ax2.set_ylabel("Recovery rate")
ax2.set_title("Fraction of trials recovered (SNR + shape veto survives near truth)")
ax2.legend(fontsize=8)
fig2.tight_layout()
fig2.savefig("injection_montecarlo_recovery_rate_2025_12_16.png", dpi=200)

# ─── Plot 3: false-alarm counts per method (uninjected series, once per radius) ──
fig3, ax3 = plt.subplots(figsize=(9, 5))
for i, method in enumerate(METHODS):
    vals = [noise_df[(noise_df.method == method) & (noise_df.radius_m == radius)]["n_false_alarms"].iloc[0] for radius in RADII_M]
    ax3.bar(x + (i - (n_methods - 1) / 2) * width, vals, width=width * 0.9, color=colors[i], label=method)
ax3.set_xticks(x)
ax3.set_xticklabels(RADII_M)
ax3.set_xlabel("KBO radius (m)")
ax3.set_ylabel("False alarms (post shape-veto, uninjected series)")
ax3.set_title("Blind-search false-alarm count per method (no injected signal present)")
ax3.legend(fontsize=8)
fig3.tight_layout()
fig3.savefig("injection_montecarlo_false_alarms_2025_12_16.png", dpi=200)

# ─── Summary table ─────────────────────────────────────────────────────────
print(f"\n{'method':18s} {'radius':>7s} {'n':>4s} {'median SNR':>11s} {'IQR':>17s} {'recovery rate':>14s}")
for method in METHODS:
    for radius in RADII_M:
        sub = trial_df[(trial_df.method == method) & (trial_df.radius_m == radius)]
        vals = sub["recovered_snr"].dropna()
        if len(vals):
            med = vals.median()
            q1, q3 = vals.quantile([0.25, 0.75])
            iqr_str = f"[{q1:6.2f}, {q3:6.2f}]"
        else:
            med, iqr_str = float("nan"), "--"
        rate = sub["recovered"].mean()
        print(f"{method:18s} {radius:7.0f} {len(sub):4d} {med:11.2f} {iqr_str:>17s} {100 * rate:13.1f}%")

print("\nSaved: injection_montecarlo_trials_2025_12_16.csv, injection_montecarlo_noise_profiles_2025_12_16.csv, "
      "injection_montecarlo_snr_boxplot_2025_12_16.png, injection_montecarlo_recovery_rate_2025_12_16.png, "
      "injection_montecarlo_false_alarms_2025_12_16.png")
