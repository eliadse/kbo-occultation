# examples/studies/synthetic_min_size.py
"""
Minimum-size pipeline on synthetic runs: injection -> blind search -> false-alarm thresholds -> efficiency -> R50.

Instead of a logged star run, the noise is generated (kbo_occultation.analysis.synthetic): Gaussian, stationary,
independent between the three channels (MAGIC-2 = A, LST-1 = B, MAGIC-1 = C), in one of two models

  white      the white level of the noise model (Instrument.noise_budget) of the star of each magnitude, in the
             conditions of a reference run, or a given SNR (--snr per --snr-dt seconds);
  coloured   the same plus the measured excess PSD of the reference run, i.e. the PSD the prediction assumes.

Every run goes through the same code as examples/studies/min_size_calibration.py (min_size.calibrate_run): the same
star-only event is injected into all channels at random times, searched blindly as the real data are, and the
detections counted for each telescope and for 2 of 3 in coincidence. The false-alarm thresholds come from the blind
search of the uninjected synthetic runs themselves (pooled over the magnitudes of a noise model and correction), not
from real_candidates.csv. The measured efficiency and R50 are compared with the prediction of analysis.min_size for
(a) the perfectly matched template and (b) the best template of the search bank, and the fitted SNR scale
measured / predicted is reported: ~1 against (b) means that the predictor is right for noise that satisfies its
assumptions (so the offsets of ``SNR_CALIBRATION`` come from the real noise), otherwise the search itself differs
(detrending, normalisation). The bank is each star's own (``mimic.star_bank``, the size of the synthetic star), as in the
real search by default, so (a) and (b) differ only by the radii of the bank (50, 100, 300, 1000 m); ``--fixed-bank`` uses the
0.03 mas bank of the older search, which loses SNR on the smaller (fainter) stars.

Cost: every (run, correction) is ~ trials x radii x 3 searches of 2 minutes of data (about 0.5 s each); per-run results
are written to <out>/trials/ and reused on a rerun.

Usage:
    python examples/studies/synthetic_min_size.py [--noise white coloured] [--mags 7 9 1] [--runs-per-mag 2]
        [--run-s 600] [--trials 30] [--radii 20 300 12] [--corrections whitened] [--fars 1 0.1 0.01]
        [--snr 10 --snr-dt 0.00524] [--profiles outputs/min_size/profiles_<nights>.npz]
        [--reference median|all|NIGHT/OBS] [--default-conditions] [--fixed-bank] [--out outputs/min_size/synthetic]
Smoke test:
    python examples/studies/synthetic_min_size.py --noise white --mags 8 8 1 --runs-per-mag 1 --run-s 300 \
        --trials 5 --radii 50 300 4 --corrections whitened --default-conditions --out outputs/min_size/synthetic_smoke
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
from kbo_occultation.analysis import synthetic as sy
from kbo_occultation.analysis.mimic import MimicConfig, observed_star_mas
from kbo_occultation.analysis.noise_removal import ComparisonConfig

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--noise", nargs="+", default=["white", "coloured"], choices=["white", "coloured"])
ap.add_argument("--mags", nargs=3, type=float, default=[7.0, 9.0, 1.0], metavar=("MIN", "MAX", "STEP"))
ap.add_argument("--runs-per-mag", type=int, default=2, help="synthetic runs per magnitude (each a new noise draw)")
ap.add_argument("--run-s", type=float, default=600.0, help="length of a synthetic run [s]")
ap.add_argument("--trials", type=int, default=30, help="injections per run, channel and radius")
ap.add_argument("--radii", nargs=3, type=float, default=[20.0, 300.0, 12], metavar=("MIN", "MAX", "N"),
                help="KBO radius grid [m], geometric")
ap.add_argument("--corrections", nargs="+", default=["whitened"], help="raw | whitened | butter_<fc>")
ap.add_argument("--fars", nargs="+", type=float, default=[1.0, 0.1, 0.01], help="coincidence FARs [1/h], first = headline")
ap.add_argument("--snr", type=float, default=None, help="photometric SNR of the star per --snr-dt (replaces the noise model's white level)")
ap.add_argument("--snr-dt", type=float, default=None, help="bin length [s] of --snr")
ap.add_argument("--profiles", default=os.path.join("outputs", "min_size", "profiles_2026-03-08_2026-08-19.npz"),
                help="noise profiles .npz of real runs: conditions of all models and the excess of 'coloured'")
ap.add_argument("--reference", default="median", help="reference run(s): median | all | NIGHT/OBSERVATION")
ap.add_argument("--default-conditions", action="store_true", help="white: do not copy the conditions of a reference run")
ap.add_argument("--fixed-bank", action="store_true", help="search with the fixed 0.03 mas bank instead of each star's own")
ap.add_argument("--out", default=os.path.join("outputs", "min_size", "synthetic"))
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()

OUT = args.out
os.makedirs(os.path.join(OUT, "trials"), exist_ok=True)
FARS = list(args.fars)
HEAD_FAR = FARS[0]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#8a5cd0"]        # one colour per noise model / night


def say(msg):
    print(msg, flush=True)


def read_csv_or_empty(path):
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def cat(frames):
    frames = [f for f in frames if len(f)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def binomial_err(hits, n):
    p = np.clip(hits / n, 0.5 / n, 1.0 - 0.5 / n)
    return np.sqrt(p * (1.0 - p) / n)


if (args.snr is None) != (args.snr_dt is None):
    sys.exit("--snr and --snr-dt go together")

cfg = ms.MinSizeConfig(time_resolved=False)               # stationary noise: one PSD per run
ccfg = ComparisonConfig(mimic=MimicConfig(per_star_size=not args.fixed_bank), n_trials=args.trials, seed=args.seed)
tol = cfg.tolerance_s()
radii = [float(r) for r in np.round(np.geomspace(args.radii[0], args.radii[1], int(args.radii[2])), 2)]
mags = [float(m) for m in np.arange(args.mags[0], args.mags[1] + 1e-9, args.mags[2])]

# ───────────────────────────── reference runs and synthetic runs ────────────
references = None
if os.path.exists(args.profiles):
    references = sy.reference_runs(ms.load_profiles(args.profiles), args.reference)
    say("reference runs: " + ", ".join(f"{g['A'].night}/{g['A'].observation}" for g in references))
elif "coloured" in args.noise:
    sys.exit(f"{args.profiles}: no noise profiles (python examples/studies/min_detectable_size.py writes them)")
else:
    say(f"{args.profiles}: not found, default conditions (airmass {cfg.airmass}, throughput 1, PSF containment 1)")

runs, mag_of = [], {}
for noise in args.noise:
    refs = references if (noise == "coloured" or (references and not args.default_conditions)) else None
    these = sy.synthetic_runs(mags, cfg, noise=noise, references=refs, snr=args.snr, snr_dt_s=args.snr_dt,
                              n_runs=args.runs_per_mag, run_s=args.run_s, seed=args.seed)
    runs += these
    mag_of.update({(n.name, o.name): o.mag_B for n, o, _ in these})
nights = sorted({n.name for n, _, _ in runs})
say(f"{len(runs)} synthetic runs ({', '.join(nights)}), B = {mags}, radii {radii[0]:g}-{radii[-1]:g} m x {len(radii)}")

# ───────────────────────────── injections (resumable) ───────────────────────
bank = None                                    # None: calibrate_run builds each star's own (mimic.star_bank)
if args.fixed_bank:
    from kbo_occultation.analysis.mimic import search_templates                   # noqa: E402
    from kbo_occultation.simulation.kinematics import EARTH_SEMIMAJOR_AXIS_AU, shadow_velocity_opposition   # noqa: E402
    bank = search_templates(ccfg.mimic, shadow_velocity_opposition(ccfg.mimic.distance_au + EARTH_SEMIMAJOR_AXIS_AU))
all_parts = {k: [] for k in ("trials", "cands", "calib", "full")}
for night, obs, group in runs:
    patterns = {}
    for corr in args.corrections:
        stem = os.path.join(OUT, "trials", f"{night.name}_{obs.name}_{corr}" + ("_fixedbank" if args.fixed_bank else ""))
        paths = {k: f"{stem}_{k}.csv" for k in all_parts}
        if all(os.path.exists(p) for p in paths.values()):
            parts = {k: read_csv_or_empty(p) for k, p in paths.items()}
        else:
            say(f"{night.name} {obs.name} {corr}: {args.trials} trials x {len(radii)} radii")
            full = []
            t, c, k = ms.calibrate_run(night, obs, group, cfg, ccfg, corr, radii, sf_mode="measured",
                                       n_trials=args.trials, seed=args.seed, bank=bank, patterns=patterns,
                                       full_candidates=full)
            parts = {"trials": t, "cands": c, "calib": k, "full": pd.DataFrame(full)}
            for key, df in parts.items():
                df.to_csv(paths[key], index=False)
        for key, df in parts.items():
            all_parts[key].append(df)
    obs.release()

trials, cands, run_calib, full = (cat(all_parts[k]) for k in ("trials", "cands", "calib", "full"))
if trials.empty:
    sys.exit("no injections (runs too short for the +-60 s windows?)")
for df in (cands, run_calib, full):
    if len(df):
        df["mag_B"] = [mag_of[(n, o)] for n, o in zip(df["night"], df["observation"])]
run_calib.to_csv(os.path.join(OUT, "run_calibration.csv"), index=False)
say("full-run search statistics (snr_std, snr_mad_sigma ~ 1 where the empirical normalisation is right):")
say(run_calib.groupby(["night", "correction", "channel"])[["snr_std", "snr_mad_sigma", "n_candidates"]]
    .mean().round(3).to_string())

# ───────────────────────────── thresholds from the synthetic search ─────────
curves = sy.rate_curves_from_candidates(full, run_calib) if len(full) else {}
thr_tables = []
for (night, corr), cv in curves.items():
    thr_tables.append(ms.thresholds_table(cv, tol, FARS).assign(night=night, correction=corr))
thresholds = cat(thr_tables)
thresholds.to_csv(os.path.join(OUT, "thresholds.csv"), index=False)
say("per-telescope SNR thresholds (2 of 3 in coincidence):")
say(thresholds[thresholds["config"].str.contains("of")].round(2).to_string(index=False))

# ───────────────────────────── measured and predicted efficiency ────────────


def predicted_with_bank():
    """Expected SNR found with the search bank of each run's star (the fixed one with --fixed-bank)."""
    banks, parts = {}, []
    for run in runs:
        obs = run[1]
        size = round(observed_star_mas(obs), 6)
        if size not in banks:
            banks[size] = ms.bank_deltas(cfg) if args.fixed_bank else \
                ms.bank_deltas(cfg, star=ms.observed_target(obs, observed_star_mas(obs)))
        parts.append(ms.predicted_snr([run], cfg, args.corrections, radii, bank=banks[size]))
    return cat(parts)


snr_pred = {"matched": ms.predicted_snr(runs, cfg, args.corrections, radii), "bank": predicted_with_bank()}
for label, df in snr_pred.items():
    df.to_csv(os.path.join(OUT, f"predicted_snr_{label}.csv"), index=False)

eff_parts = []
for (night, corr), g in trials.groupby(["night", "correction"]):
    thr = thresholds[(thresholds["night"] == night) & (thresholds["correction"] == corr)]
    sel_c = cands[(cands["night"] == night) & (cands["correction"] == corr)] if len(cands) else cands
    meas = ms.measured_efficiency(g, sel_c, thr, tol)
    for label, df in snr_pred.items():
        sp = df[(df["night"] == night) & (df["correction"] == corr)]
        pred = ms.predicted_efficiency(sp, thr).rename(columns={"eff_pred": f"eff_pred_{label}"})
        meas = meas.merge(pred[["night", "observation", "correction", "radius_m", "config", "far_per_hr",
                                f"eff_pred_{label}"]],
                          on=["night", "observation", "correction", "radius_m", "config", "far_per_hr"], how="left")
    eff_parts.append(meas)
eff = cat(eff_parts)
eff["mag_B"] = [mag_of[(n, o)] for n, o in zip(eff["night"], eff["observation"])]
eff.to_csv(os.path.join(OUT, "efficiency_per_run.csv"), index=False)

# pooled over the runs of a magnitude: hits / n, the mean prediction
key = ["night", "correction", "mag_B", "radius_m", "config", "far_per_hr"]
pooled = (eff.groupby(key).agg(hits=("hits", "sum"), n=("n", "sum"), threshold=("threshold", "first"),
                               eff_pred_matched=("eff_pred_matched", "mean"), eff_pred_bank=("eff_pred_bank", "mean"))
          .reset_index())
pooled["eff"] = pooled["hits"] / pooled["n"]
pooled["eff_err"] = binomial_err(pooled["hits"].to_numpy(float), pooled["n"].to_numpy(float))
pooled.to_csv(os.path.join(OUT, "efficiency.csv"), index=False)

by = ["night", "correction", "mag_B", "config", "far_per_hr"]
r50 = None
for col, name in (("eff", "meas"), ("eff_pred_matched", "matched"), ("eff_pred_bank", "bank")):
    part = ms.efficiency_r50(pooled, col, by).rename(columns={"r50_m": f"r50_{name}", "status": f"status_{name}"})
    r50 = part if r50 is None else r50.merge(part, on=by)
r50["ratio_bank"] = r50["r50_bank"] / r50["r50_meas"]
r50.to_csv(os.path.join(OUT, "r50.csv"), index=False)
say(f"R50 [m] at FAR = {HEAD_FAR:g}/h (measured / matched template / search bank):")
say(r50[(r50["far_per_hr"] == HEAD_FAR)].sort_values(["night", "correction", "config", "mag_B"])
    [["night", "correction", "config", "mag_B", "r50_meas", "r50_matched", "r50_bank", "status_meas", "ratio_bank"]]
    .round(2).to_string(index=False))

# ───────────────────────────── SNR scale: measured / predicted ──────────────
rows, dev = [], []
for label, sp in snr_pred.items():
    per_run, _ = ms.fit_snr_calibration(trials, sp, run_calib)
    per_run.assign(prediction=label).to_csv(os.path.join(OUT, f"snr_calibration_per_run_{label}.csv"), index=False)
    s = (per_run.groupby(["night", "correction", "kind"])
         .agg(scale=("scale", "median"), scale_min=("scale", "min"), scale_max=("scale", "max"),
              scatter=("scatter", "median"), n_run_channels=("scale", "size")).reset_index())
    rows.append(s.assign(prediction=label))
    for (night, corr), g in eff.groupby(["night", "correction"]):
        thr = thresholds[(thresholds["night"] == night) & (thresholds["correction"] == corr)]
        sps = sp[(sp["night"] == night) & (sp["correction"] == corr)]
        try:
            best, _ = ms.scale_by_deviance(g, sps, thr, corr, scatter=1.0)
        except (ValueError, KeyError, IndexError):
            best = np.nan
        dev.append({"night": night, "correction": corr, "prediction": label, "scale_deviance": best})
calibration = cat(rows).merge(pd.DataFrame(dev), on=["night", "correction", "prediction"], how="left")
calibration.to_csv(os.path.join(OUT, "snr_calibration.csv"), index=False)
say("SNR scale (measured / predicted peak SNR, median over run-channels; scale_deviance = likelihood fit on the "
    "2-of-3 detections):")
say(calibration.round(3).to_string(index=False))
say("for comparison, the real runs: SNR_CALIBRATION_STATIONARY = " + repr(ms.SNR_CALIBRATION_STATIONARY))

# ───────────────────────────── figures ──────────────────────────────────────
color_of = {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(nights)}
for corr in args.corrections:
    sub = pooled[(pooled["correction"] == corr) & (pooled["config"] == "2of3") & (pooled["far_per_hr"] == HEAD_FAR)]
    if sub.empty:
        continue
    fig, axes = plt.subplots(len(nights), len(mags), figsize=(3.2 * len(mags) + 0.6, 2.5 * len(nights) + 0.9),
                             sharex=True, sharey=True, squeeze=False)
    for i, night in enumerate(nights):
        for j, mag in enumerate(mags):
            ax = axes[i, j]
            style(ax)
            g = sub[(sub["night"] == night) & (sub["mag_B"] == mag)].sort_values("radius_m")
            ax.errorbar(g["radius_m"], g["eff"], g["eff_err"], fmt="o", ms=4, color=color_of[night], elinewidth=1,
                        label="injections")
            ax.plot(g["radius_m"], g["eff_pred_matched"], "--", color=INK, lw=1.2, label="matched template")
            ax.plot(g["radius_m"], g["eff_pred_bank"], "-", color=INK, lw=1.5, label="search bank")
            ax.axhline(0.5, color=GRID, lw=1)
            ax.set_xscale("log")
            ax.set_ylim(-0.05, 1.05)
            ticks = [t for t in (20, 50, 100, 200) if radii[0] / 1.2 <= t <= radii[-1] * 1.2]
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(t) for t in ticks])
            ax.minorticks_off()
            if i == 0:
                ax.set_title(f"B = {mag:g}", fontsize=9, color=INK)
            if j == 0:
                ax.set_ylabel(f"{night}\nefficiency", fontsize=8, color=MUTED)
            if i == len(nights) - 1:
                ax.set_xlabel("KBO radius [m]", fontsize=8, color=MUTED)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(f"2 of 3 in coincidence at FAR = {HEAD_FAR:g}/h on synthetic runs, correction: {corr}",
                 fontsize=10, color=INK)
    fig.tight_layout(rect=[0, 0.04, 1, 0.97])
    fig.savefig(os.path.join(OUT, f"efficiency_{corr}.png"), dpi=140)
    plt.close(fig)

for corr in args.corrections:
    sub = r50[(r50["correction"] == corr) & (r50["config"] == "2of3") & (r50["far_per_hr"] == HEAD_FAR)]
    if sub.empty:
        continue
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    style(ax)
    for night in nights:
        g = sub[sub["night"] == night].sort_values("mag_B")
        ax.plot(g["mag_B"], g["r50_meas"], "o", ms=5, color=color_of[night], label=f"{night}: injections")
        ax.plot(g["mag_B"], g["r50_matched"], "--", color=color_of[night], lw=1.2)
        ax.plot(g["mag_B"], g["r50_bank"], "-", color=color_of[night], lw=1.5)
    ax.plot([], [], "--", color=INK, lw=1.2, label="matched template")
    ax.plot([], [], "-", color=INK, lw=1.5, label="search bank")
    ax.set_xlabel("B magnitude", fontsize=8, color=MUTED)
    ax.set_ylabel("R50 [m]", fontsize=8, color=MUTED)
    ax.set_title(f"R50, 2 of 3, FAR = {HEAD_FAR:g}/h, correction: {corr}", fontsize=9, color=INK)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"r50_vs_mag_{corr}.png"), dpi=140)
    plt.close(fig)
say(f"wrote {OUT}/")
