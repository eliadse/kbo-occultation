# examples/studies/min_size_calibration.py
"""
Calibrate the minimum-size predictor (kbo_occultation.analysis.min_size) with array injections.

For every logged star run with >= 2 channels the same star-only event is injected into all channels
(shadow offsets between the telescopes), corrected, blind-searched as the real data are, and the
detections are counted for each single telescope and for k of n in coincidence at the real-data
thresholds. The measured efficiency and peak SNR are compared with the step-1 prediction; the fitted
SNR scale / scatter go to ``snr_calibration.csv`` (copy them into ``min_size.SNR_CALIBRATION``).

Needs ``outputs/noise_removal/real_candidates.csv`` (examples/studies/noise_removal_comparison.py) with
the nights and corrections used here. Per-run results are written to ``<out>/trials/`` as they finish
and reused on a rerun, so an interrupted run resumes.

Usage:
    python examples/studies/min_size_calibration.py [NIGHT ...] [--trials 50]
        [--corrections whitened raw butter_1] [--radii 50 100 300] [--model-sf-trials 20]
Smoke test:
    python examples/studies/min_size_calibration.py 2026-08-19 --trials 5 --corrections whitened \
        --radii 50 100 300 --model-sf-trials 0 --shifts 200 --out outputs/min_size/calibration_smoke
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from kbo_occultation.analysis import min_size as ms
from kbo_occultation.analysis.noise_removal import ComparisonConfig
from kbo_occultation.config import ORM_ARRAY
from kbo_occultation.observations.log import load_observation_log

ap = argparse.ArgumentParser()
ap.add_argument("nights", nargs="*", default=["2026-03-08", "2026-08-19"])
ap.add_argument("--trials", type=int, default=50, help="injections per run, channel and radius")
ap.add_argument("--corrections", nargs="+", default=["whitened", "raw", "butter_1"])
ap.add_argument("--radii", nargs="+", type=float, default=list(ms.CALIBRATION_RADII_M))
ap.add_argument("--model-sf-trials", type=int, default=20,
                help="trials of the model-star-fraction variant (whitened only; 0 = skip)")
ap.add_argument("--shifts", type=int, default=2000, help="time slides for the coincidence FAR check")
ap.add_argument("--candidates", default=os.path.join("outputs", "noise_removal", "real_candidates.csv"))
ap.add_argument("--profiles", default=None, help="noise profiles .npz (default: outputs/min_size/profiles_<nights>.npz)")
ap.add_argument("--out", default=os.path.join("outputs", "min_size", "calibration"))
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--include-excluded", action="store_true")
ap.add_argument("--force", action="store_true", help="go on although the candidate list looks stale")
args = ap.parse_args()

OUT = args.out
os.makedirs(os.path.join(OUT, "trials"), exist_ok=True)
FARS = [1.0, 0.1, 0.01]
COLORS = {"MAGIC-2": "#2a78d6", "MAGIC-1": "#eb6834", "LST-1": "#1baf7a"}     # validated categorical slots 1-3
MARKERS = {"MAGIC-2": "o", "MAGIC-1": "s", "LST-1": "^"}                       # secondary encoding (aqua < 3:1)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def say(msg):
    print(msg, flush=True)


def read_csv_or_empty(path):
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


log = load_observation_log()
cfg = ms.MinSizeConfig()
ccfg = ComparisonConfig(n_trials=args.trials, seed=args.seed)
tol = cfg.tolerance_s()
radii = [float(r) for r in args.radii]

# ───────────────────────────── noise profiles and runs ──────────────────────
prof_path = args.profiles or os.path.join("outputs", "min_size", "profiles_" + "_".join(args.nights) + ".npz")
if os.path.exists(prof_path):
    profiles = ms.load_profiles(prof_path)
    say(f"profiles: {len(profiles)} from {prof_path}")
else:
    say("measuring noise profiles (slow: loads the light curves)")
    profiles = ms.measure_noise_profiles(log, args.nights, cfg, include_excluded=args.include_excluded, progress=say)
    os.makedirs(os.path.dirname(prof_path), exist_ok=True)
    ms.save_profiles(profiles, prof_path)
runs = ms.calibration_runs(log, args.nights, profiles, include_excluded=args.include_excluded)
if not runs:
    sys.exit("no run with >= 2 channels and noise profiles in the selected nights")
say("runs: " + ", ".join(f"{n.name}/{o.name}({len(g)}ch)" for n, o, g in runs))
mag = {(n.name, o.name): o.mag_B for n, o, _ in runs}

# ───────────────────────────── candidate list check, thresholds, FAR ────────
cand_df = pd.read_csv(args.candidates, usecols=["night", "observation", "method"])
problems = []
seen_runs = set(zip(cand_df["night"], cand_df["observation"]))          # a run with no candidates under a method is fine
for n, o, _ in runs:
    if (n.name, o.name) not in seen_runs:
        problems.append(f"{n.name}/{o.name}: no candidates under any method (night or run missing)")
for corr in args.corrections:
    name = ms.csv_method_name(corr)
    nights_with = set(cand_df.loc[cand_df["method"] == name, "night"])
    problems += [f"{name}: no candidates in {n}" for n in args.nights if n not in nights_with]
thr_csv = os.path.join(os.path.dirname(args.candidates), "thresholds.csv")
if os.path.exists(thr_csv):
    th = pd.read_csv(thr_csv).set_index("method")
    for corr in args.corrections:
        name = ms.csv_method_name(corr)
        curves_all = ms.rate_curves(args.candidates, log, name, None, include_excluded=args.include_excluded)
        hours = sum(c.hours for c in curves_all.values())
        if name in th.index and abs(hours - th.loc[name, "hours"]) > 0.01 * th.loc[name, "hours"]:
            problems.append(f"{name}: {hours:.3f} telescope-h in the log vs {th.loc[name, 'hours']:.3f} in thresholds.csv")
if problems:
    say("candidate list looks stale:\n  " + "\n  ".join(problems) +
        "\nregenerate it: python examples/studies/noise_removal_comparison.py --methods raw whitened butter_1Hz")
    if not args.force:
        sys.exit(1)

curves, far_checks, thr_tables = {}, {}, []
for corr in args.corrections:
    name = ms.csv_method_name(corr)
    curves[corr] = ms.rate_curves(args.candidates, log, name, args.nights, include_excluded=args.include_excluded)
    far_checks[corr] = ms.timeslide_far(args.candidates, log, name, args.nights, tol, n_shifts=args.shifts,
                                        curves=curves[corr], seed=args.seed, include_excluded=args.include_excluded)
    far_checks[corr].assign(correction=corr).to_csv(os.path.join(OUT, f"far_check_{corr}.csv"), index=False)
    thr_tables.append(ms.thresholds_table(curves[corr], tol, FARS, far_checks[corr]).assign(correction=corr))
thresholds = pd.concat(thr_tables, ignore_index=True)
thresholds.to_csv(os.path.join(OUT, "thresholds.csv"), index=False)
say(thresholds[thresholds["config"].str.contains("of")].round(2).to_string(index=False))

# ───────────────────────────── injections (resumable) ───────────────────────
bank = None
all_trials, all_cands, all_calib = [], [], []


def run_injections(corr, sf_mode, n_trials, night, obs, group, patterns):
    global bank
    stem = os.path.join(OUT, "trials", f"{night.name}_{obs.name}_{corr}_{sf_mode}")
    paths = [f"{stem}_{k}.csv" for k in ("trials", "cands", "calib")]
    if all(os.path.exists(p) for p in paths):
        return [read_csv_or_empty(p) for p in paths]
    if bank is None:
        from kbo_occultation.analysis.mimic import search_templates
        from kbo_occultation.simulation.kinematics import EARTH_SEMIMAJOR_AXIS_AU, shadow_velocity_opposition
        bank = search_templates(ccfg.mimic, shadow_velocity_opposition(ccfg.mimic.distance_au + EARTH_SEMIMAJOR_AXIS_AU))
    say(f"{night.name} {obs.name} {corr} [{sf_mode}]: {n_trials} trials x {len(radii)} radii")
    out = ms.calibrate_run(night, obs, group, cfg, ccfg, corr, radii, sf_mode=sf_mode, n_trials=n_trials,
                           seed=args.seed, bank=bank, patterns=patterns)
    for df, p in zip(out, paths):
        df.to_csv(p, index=False)
    return list(out)


for night, obs, group in runs:
    patterns = {}
    for corr in args.corrections:
        for df, store in zip(run_injections(corr, "measured", args.trials, night, obs, group, patterns),
                             (all_trials, all_cands, all_calib)):
            store.append(df)
    if args.model_sf_trials > 0 and "whitened" in args.corrections:
        for df, store in zip(run_injections("whitened", "model", args.model_sf_trials, night, obs, group, patterns),
                             (all_trials, all_cands, all_calib)):
            store.append(df)

trials, cands, run_calib = (pd.concat([d for d in s if len(d)], ignore_index=True) if any(len(d) for d in s)
                            else pd.DataFrame() for s in (all_trials, all_cands, all_calib))
if trials.empty:
    sys.exit("no injections (good-time windows too short?)")
run_calib.to_csv(os.path.join(OUT, "run_calibration.csv"), index=False)
real_check = run_calib[run_calib["sf_mode"] == "measured"].drop_duplicates(["night", "observation", "correction", "channel"])
say("full-run search statistics (snr_std ~ factor the robust normalisation gives over an rms prediction):")
say(real_check.round(3).to_string(index=False))

# ───────────────────────────── efficiency: measured vs predicted ────────────
snr_pred = ms.predicted_snr(runs, cfg, args.corrections, radii)
snr_pred.to_csv(os.path.join(OUT, "predicted_snr.csv"), index=False)
meas_parts = []
for corr in args.corrections:
    t, c = trials[trials["correction"] == corr], cands[cands["correction"] == corr] if len(cands) else cands
    meas_parts.append(ms.measured_efficiency(t, c, thresholds[thresholds["correction"] == corr], tol))
eff_meas = pd.concat(meas_parts, ignore_index=True)

per_run, summary = ms.fit_snr_calibration(trials, snr_pred, run_calib)
per_run.to_csv(os.path.join(OUT, "snr_calibration_per_run.csv"), index=False)
scales = {}
dev = []
for corr in args.corrections:
    row = summary[(summary["correction"] == corr) & (summary["sf_mode"] == "measured") & (summary["kind"] == "all")]
    if row.empty:
        continue
    scale, scatter = float(row["scale"].iloc[0]), float(row["scatter"].iloc[0])
    scatter = scatter if np.isfinite(scatter) and scatter > 0 else 1.0
    scales[corr] = (scale, scatter)
    m = eff_meas[eff_meas["sf_mode"] == "measured"]
    best, _ = ms.scale_by_deviance(m, snr_pred, thresholds[thresholds["correction"] == corr], corr, scatter=scatter)
    dev.append({"correction": corr, "sf_mode": "measured", "kind": "all", "scale_deviance": best})
summary = summary.merge(pd.DataFrame(dev), on=["correction", "sf_mode", "kind"], how="left") if dev else summary
summary.to_csv(os.path.join(OUT, "snr_calibration.csv"), index=False)
say("SNR calibration (measured / predicted), median over runs:")
say(summary.round(3).to_string(index=False))
say("SNR_CALIBRATION = " + repr({k: (round(a, 3), round(b, 3)) for k, (a, b) in scales.items()}))

keys = ["night", "observation", "correction", "radius_m", "config", "far_per_hr"]
pred0 = pd.concat([ms.predicted_efficiency(snr_pred[snr_pred["correction"] == c],
                                           thresholds[thresholds["correction"] == c]) for c in args.corrections])
pred1 = pd.concat([ms.predicted_efficiency(snr_pred[snr_pred["correction"] == c],
                                           thresholds[thresholds["correction"] == c], scales) for c in args.corrections])
eff = (eff_meas.merge(pred0[keys + ["eff_pred"]], on=keys, how="left")
       .merge(pred1[keys + ["eff_pred"]].rename(columns={"eff_pred": "eff_pred_cal"}), on=keys, how="left"))
eff.to_csv(os.path.join(OUT, "efficiency.csv"), index=False)

# ───────────────────────────── R50 residuals per run ────────────────────────
e1 = eff[(eff["sf_mode"] == "measured")]
by = ["night", "observation", "correction", "config", "far_per_hr"]
r_meas = ms.efficiency_r50(e1, "eff", by).rename(columns={"r50_m": "r50_meas", "status": "status_meas"})
r_pred = ms.efficiency_r50(e1, "eff_pred", by).rename(columns={"r50_m": "r50_pred", "status": "status_pred"})
r_cal = ms.efficiency_r50(e1, "eff_pred_cal", by).rename(columns={"r50_m": "r50_cal", "status": "status_cal"})
r50 = r_meas.merge(r_pred, on=by).merge(r_cal, on=by)
r50["resid_pred"] = r50["r50_pred"] / r50["r50_meas"] - 1.0
r50["resid_cal"] = r50["r50_cal"] / r50["r50_meas"] - 1.0
r50.to_csv(os.path.join(OUT, "r50_residuals.csv"), index=False)
say("R50 [m] per run, 2 of n at FAR 1/h (measured / step-1 / calibrated):")
say(r50[r50["config"].str.contains("of") & (r50["far_per_hr"] == 1.0)]
    [["night", "observation", "correction", "config", "r50_meas", "r50_pred", "r50_cal", "resid_pred", "resid_cal"]]
    .round(2).to_string(index=False))

ex = ms.excess_scaling_test(profiles)
ex.to_csv(os.path.join(OUT, "excess_scaling.csv"), index=False)
say("excess noise vs star fraction (slope 0: fixed in star units, -2: fixed in total-signal units):")
say(ex.round(3).to_string(index=False))

# ───────────────────────────── figures ──────────────────────────────────────
for corr in args.corrections:
    sub = eff[(eff["correction"] == corr) & (eff["sf_mode"] == "measured") & (eff["far_per_hr"] == 1.0)]
    run_keys = sorted(set(zip(sub["night"], sub["observation"])))
    if not run_keys:
        continue
    cfgs = [ORM_ARRAY.by_channel(c).name for c in "ACB"] + ["2of3"]
    fig, axes = plt.subplots(len(run_keys), 4, figsize=(13, 2.3 * len(run_keys) + 0.8), sharex=True, sharey=True,
                             squeeze=False)
    for i, (night, obs) in enumerate(run_keys):
        for j, label in enumerate(cfgs):
            ax = axes[i, j]
            style(ax)
            config = "2of3" if label == "2of3" else next(c for c in "ABC" if ORM_ARRAY.by_channel(c).name == label)
            g = sub[(sub["night"] == night) & (sub["observation"] == obs) & (sub["config"] == config)].sort_values("radius_m")
            color = INK if label == "2of3" else COLORS[label]
            ax.errorbar(g["radius_m"], g["eff"], g["eff_err"], fmt=MARKERS.get(label, "D"), ms=4, color=color,
                        elinewidth=1, label="injections")
            ax.plot(g["radius_m"], g["eff_pred"], "--", color=color, lw=1.2, label="step 1 (snr_scale 1)")
            ax.plot(g["radius_m"], g["eff_pred_cal"], "-", color=color, lw=1.5, label="calibrated")
            ax.set_xscale("log")
            ax.set_xlim(min(radii) / 1.2, max(radii) * 1.2)
            ticks = [t for t in (50, 100, 200, 400) if min(radii) / 1.2 <= t <= max(radii) * 1.2]
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(t) for t in ticks])
            ax.minorticks_off()
            ax.set_ylim(-0.05, 1.05)
            if i == 0:
                ax.set_title(label if label != "2of3" else "2 of 3 in coincidence", fontsize=9, color=INK)
            if j == 0:
                ax.set_ylabel(f"{obs}\nB={mag[(night, obs)]:.2f}\nefficiency", fontsize=8, color=MUTED)
            if i == len(run_keys) - 1:
                ax.set_xlabel("KBO radius [m]", fontsize=8, color=MUTED)
    handles, labels = axes[0, 3].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(f"Detection efficiency at FAR = 1/h, correction: {corr}", fontsize=10, color=INK)
    fig.tight_layout(rect=[0, 0.025, 1, 0.975])
    fig.savefig(os.path.join(OUT, f"efficiency_{corr}.png"), dpi=140)
    plt.close(fig)

pr = per_run[per_run["sf_mode"] == "measured"].dropna(subset=["snr_std"]) if "snr_std" in per_run else pd.DataFrame()
if len(pr):
    corrs = [c for c in args.corrections if c in set(pr["correction"])]
    lim = [0.8, max(3.0, float(pr["snr_std"].max()) * 1.1, float(pr["scale"].max()) * 1.1)]
    fig, axes = plt.subplots(1, len(corrs), figsize=(4.4 * len(corrs), 4.2), sharex=True, sharey=True, squeeze=False)
    for ax, corr in zip(axes[0], corrs):
        style(ax)
        for tel, g in pr[pr["correction"] == corr].groupby("telescope"):
            ax.scatter(g["snr_std"], g["scale"], s=26, color=COLORS.get(tel, INK), marker=MARKERS.get(tel, "D"),
                       label=tel, edgecolor="#fcfcfb", linewidth=0.8)
        ax.plot(lim, lim, color=MUTED, lw=1, ls=":", label="scale = snr_std")
        ax.set_title(corr, fontsize=9, color=INK)
        ax.set_xlabel("SNR-trace std of the full run", fontsize=8, color=MUTED)
    axes[0, 0].set_ylabel("fitted measured / predicted SNR", fontsize=8, color=MUTED)
    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "scale_vs_snr_std.png"), dpi=140)
    plt.close(fig)

fig, axes = plt.subplots(1, len(args.corrections), figsize=(4.6 * len(args.corrections), 3.8), squeeze=False)
for ax, corr in zip(axes[0], args.corrections):
    style(ax)
    fc = far_checks[corr]
    ax.plot(fc["snr"], fc["analytic_rate_per_hr"], color=INK, lw=1.5, label="independent streams")
    ok = fc["slide_count"] > 0
    ax.errorbar(fc["snr"][ok], fc["slide_rate_per_hr"][ok], fc["slide_rate_per_hr"][ok] / np.sqrt(fc["slide_count"][ok]),
                fmt="o", ms=3.5, color=COLORS["MAGIC-2"], elinewidth=1, label="time slides")
    z = fc["zero_lag_count"] > 0
    ax.plot(fc["snr"][z], fc["zero_lag_rate_per_hr"][z], "s", ms=4, color=COLORS["MAGIC-1"], label="zero lag")
    for far in FARS:
        ax.axhline(far, color=GRID, lw=1)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1e3)
    ax.set_xlabel("weakest-member SNR", fontsize=8, color=MUTED)
    ax.set_ylabel("coincidence rate [1/h]", fontsize=8, color=MUTED)
    ax.set_title(corr, fontsize=9, color=INK)
axes[0, 0].legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "coincidence_far.png"), dpi=140)
plt.close(fig)
say(f"wrote {OUT}/")
