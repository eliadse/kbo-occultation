# examples/studies/catalogue_min_size.py
"""
Smallest detectable KBO around the top-scored stars of the observing catalogue (kbo_occultation.analysis.min_size).

The catalogue (Gaia DR3 + Simbad + a score for the chance of occultation and detection) gives every star its own angular
radius and temperature, so the template of the prediction (the star seen through the telescope's band) and the colour of
the photon budget are those of the star. The ``--top`` stars by ``--sort-by`` (minus ``--exclude``) are run as such. A set
of dummy stars follows the size-magnitude-temperature relation fitted to them (``min_size.fit_star_relation``: fixed
surface brightness, theta(B) = theta_7 10^(-0.2 (B - 7)), median temperature) on a grid of magnitudes, to see R50 as a
function of the magnitude and of the size beyond the few stars that exist.

Everything else is as ``min_detectable_size.py``: R50 is the radius with 50 % detection by >= 2 of 3 telescopes at a
coincidence false-alarm rate of 1 per hour, the noise is the measured one of the logged runs (calibrated, time resolved),
"theory" is the ideal system. The star-specific template is assumed known exactly; ``bank_check.csv`` gives what the
4-template bank of the real search loses against it at R50.

Usage:
    python examples/studies/catalogue_min_size.py [NIGHT ...] [--top 5] [--exclude "HD 165132"] [--mags 7 9.5 0.25]
        [--fars 1 0.1 0.01] [--excess-scalings star total] [--out outputs/min_size/catalogue]
Smoke test:
    python examples/studies/catalogue_min_size.py --mags 7 9 1 --radii 10 1000 15 --fars 1 --excess-scalings star \
        --out outputs/min_size/catalogue_smoke
"""
import argparse
import json
import os
import sys
from dataclasses import replace

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from kbo_occultation.analysis import min_size as ms
from kbo_occultation.observations.log import load_observation_log
from kbo_occultation.simulation.instruments import Instrument
from kbo_occultation.simulation.physics import AU_m, mas_to_rad

CATALOGUE = os.path.join(os.path.dirname(ms.__file__), "..", "data", "Ecliptic_stars_noOutliers_brighter.csv")
COLUMNS = ("name=MAIN_ID,mag_B=B,angular_radius_mas=ang_radius,temperature_K=teff_gspphot,spectral_type=SP_TYPE")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("nights", nargs="*", default=["2025-12-16", "2026-03-08", "2026-08-19"],
                help="nights of the noise profiles and rate curves (the array statistics use runs with --min-channels)")
ap.add_argument("--catalogue", default=CATALOGUE)
ap.add_argument("--columns", default=COLUMNS, help="catalogue column map: field=column,...")
ap.add_argument("--sort-by", default="score", help="catalogue column, descending")
ap.add_argument("--top", type=int, default=5)
ap.add_argument("--exclude", nargs="*", default=["HD 165132"],
                help="stars left out (HD 165132: Simbad B 9.11 vs Gaia G 8.04, inconsistent with its size)")
ap.add_argument("--min-channels", type=int, default=3, help="the headline uses the runs with this many channels")
ap.add_argument("--correction", default="whitened", help="raw | whitened | butter_<fc>")
ap.add_argument("--fars", nargs="+", type=float, default=[1.0, 0.1, 0.01], help="coincidence FARs [1/h], the first is the headline")
ap.add_argument("--mags", nargs=3, type=float, default=[7.0, 9.5, 0.25], metavar=("MIN", "MAX", "STEP"),
                help="magnitudes of the dummy stars")
ap.add_argument("--radii", nargs=3, type=float, default=[10.0, 1000.0, 41], metavar=("MIN", "MAX", "N"),
                help="KBO radius grid [m], geometric")
ap.add_argument("--focus", nargs="+", type=float, default=[1.0], help="LST-1 star-light gains (2 = infinity focus)")
ap.add_argument("--excess-scalings", nargs="+", default=["star", "total"])
ap.add_argument("--stationary", action="store_true", help="one PSD per run instead of one per time block")
ap.add_argument("--profiles", default=None, help="noise profiles .npz (default outputs/min_size/profiles_<nights>.npz)")
ap.add_argument("--candidates", default=os.path.join("outputs", "noise_removal", "real_candidates.csv"))
ap.add_argument("--calibration", default=os.path.join("outputs", "min_size", "calibration"),
                help="directory of min_size_calibration.py (snr_calibration_per_run.csv gives the band)")
ap.add_argument("--generic", default=os.path.join("outputs", "min_size", "r50_summary.csv"),
                help="r50_summary.csv of min_detectable_size.py, overlaid as the generic 0.03 mas star")
ap.add_argument("--out", default=os.path.join("outputs", "min_size", "catalogue"))
ap.add_argument("--include-excluded", action="store_true")
args = ap.parse_args()

OUT = args.out
os.makedirs(OUT, exist_ok=True)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
COLORS = {"MAGIC-2": BLUE, "MAGIC-1": ORANGE, "LST-1": GREEN, "2of3": INK, "3of3": MUTED}


def say(msg):
    print(msg, flush=True)


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


log = load_observation_log()
radii = np.geomspace(args.radii[0], args.radii[1], int(args.radii[2]))
cfg = ms.MinSizeConfig(radii_m=radii, correction=args.correction, coinc_far_per_hr=args.fars[0],
                       extra_fars=tuple(args.fars[1:]), time_resolved=not args.stationary).calibrated()
say(f"correction {cfg.correction}, snr_scale {cfg.snr_scale:.3f}, scatter {cfg.snr_scatter:.3f}, "
    f"{'time resolved' if cfg.time_resolved else 'stationary'} noise, FAR {args.fars}")
HEAD_FAR = cfg.coinc_far_per_hr
tol = cfg.tolerance_s()

# ───────────────────────────── noise, thresholds ────────────────────────────
prof_path = args.profiles or os.path.join("outputs", "min_size", "profiles_" + "_".join(args.nights) + ".npz")
profiles = ms.load_or_measure_profiles(log, args.nights, cfg, prof_path, include_excluded=args.include_excluded,
                                       progress=say)
runs = ms.calibration_runs(log, args.nights, profiles, include_excluded=args.include_excluded,
                           min_channels=args.min_channels)
if not runs:
    sys.exit(f"no run with >= {args.min_channels} channels in {args.nights}")
array_profiles = [p for _, _, group in runs for p in group]
say(f"{len(runs)} runs: " + ", ".join(f"{n.name}/{o.name}" for n, o, _ in runs))

curves = ms.rate_curves(args.candidates, log, ms.csv_method_name(cfg.correction), args.nights,
                        include_excluded=args.include_excluded)
ms.thresholds_table(curves, tol, args.fars).to_csv(os.path.join(OUT, "thresholds.csv"), index=False)

band = None
per_run_csv = os.path.join(args.calibration, "snr_calibration_per_run.csv")
if os.path.exists(per_run_csv):
    try:
        band = ms.scale_quantiles(per_run_csv, cfg.correction)
        say(f"calibration band (16-84 % of the per-run scale): {band[0]:.2f} - {band[1]:.2f} (used: {cfg.snr_scale:.2f})")
    except ValueError as e:
        say(f"no calibration band: {e}")

# ───────────────────────────── targets ──────────────────────────────────────
real = ms.load_target_catalogue(args.catalogue, ms.parse_column_map(args.columns), sort_by=args.sort_by, top=args.top,
                                exclude=args.exclude)
theta7, temp = ms.fit_star_relation(real)
dummy = ms.generic_targets(np.arange(args.mags[0], args.mags[1] + 1e-9, args.mags[2]), ref_radius_mas=theta7,
                           temperature_K=temp, prefix="top B")
say(f"\ntop {len(real)} by {args.sort_by} (excluding {', '.join(args.exclude) or 'none'}):")
for t in real:
    say(f"  {t.name:12s} {t.spectral_type or '':14s} B={t.mag_B:5.2f}  theta={t.angular_radius_mas:.4f} mas  Teff={t.temperature_K:7.0f} K")
say(f"relation: theta(B) = {theta7:.4f} mas * 10^(-0.2 (B - 7)), T = {temp:.0f} K -> {len(dummy)} dummy stars")

cat = pd.read_csv(args.catalogue)
cat["name"] = [" ".join(str(n).split()) for n in cat[ms.parse_column_map(args.columns)["name"]]]
keep = ["name", "SP_TYPE", "B", "V", "phot_g_mean_mag", "teff_gspphot", "ang_radius", "radius_combined",
        "distance_gspphot", args.sort_by]
cat.set_index("name").loc[[t.name for t in real], [c for c in keep[1:] if c in cat]].reset_index() \
    .to_csv(os.path.join(OUT, "top_stars.csv"), index=False)
with open(os.path.join(OUT, "relation.json"), "w") as f:
    json.dump({"ref_mag": 7.0, "ref_radius_mas": theta7, "temperature_K": temp, "stars": [t.name for t in real],
               "excluded": args.exclude, "sort_by": args.sort_by}, f, indent=2)

# ───────────────────────────── R50 ──────────────────────────────────────────
tables, effs, summaries = {}, {}, {}
for label, stars in (("top", real), ("dummy", dummy)):
    say(f"\nR50 table, {label} stars ...")
    tables[label], effs[label] = ms.min_size_table(stars, array_profiles, cfg, curves, focus_gains=args.focus,
                                                   excess_scalings=args.excess_scalings, scale_band=band,
                                                   return_efficiency=True, progress=say)
    tables[label].to_csv(os.path.join(OUT, f"{label}_r50_table.csv"), index=False)
    summaries[label] = ms.summarize(tables[label])
    summaries[label].to_csv(os.path.join(OUT, f"{label}_r50_summary.csv"), index=False)
eff = pd.concat([e.assign(set=k) for k, e in effs.items()], ignore_index=True)
eff = eff[(eff["far_per_hr"] == HEAD_FAR) & (eff["lst_focus"] == 1.0)]
eff.groupby(["set", "target", "mag_B", "noise", "excess_scaling", "config", "radius_m"], as_index=False)["eff"] \
    .median().to_csv(os.path.join(OUT, "efficiency.csv"), index=False)

main_scaling = args.excess_scalings[0]
tel_names = sorted({p.telescope for p in array_profiles})
configs_present = tel_names + [c for c in ("2of3", "3of3") if c in set(tables["top"]["config"])]


def pick(df, noise="measured", focus=1.0, far=None, scaling=None):
    q = df[(df["noise"] == noise) & (df["lst_focus"] == focus) & (df["far_per_hr"] == (HEAD_FAR if far is None else far))]
    if noise == "measured":
        q = q[q["excess_scaling"] == (scaling or main_scaling)]
    return q


def headline(summary, order):
    s = pick(summary)
    lines = ["| star | B | theta [mas] | " + " | ".join(configs_present) + " |", "|---|---|---|" + "---|" * len(configs_present)]
    for name in order:
        cells = []
        for c in configs_present:
            r = s[(s["target"] == name) & (s["config"] == c)]
            if r.empty:
                cells.append("-")
                continue
            r = r.iloc[0]
            text = f"{r.r50_median:.0f} ({r.r50_min:.0f}-{r.r50_max:.0f})"
            lo, hi = r.get("r50_lo_median", np.nan), r.get("r50_hi_median", np.nan)
            if np.isfinite(lo) and np.isfinite(hi):
                text += f" [{lo:.0f}-{hi:.0f}]"
            cells.append(text + ("*" if r.thr_at_floor else ""))
        first = s[s["target"] == name].iloc[0]
        lines.append(f"| {name} | {first.mag_B:.2f} | {first.star_mas:.4f} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


head = (f"R50 [m], measured noise, excess '{main_scaling}', FAR {HEAD_FAR:g}/h: median (run range) [calibration band]; "
        "* = threshold at the candidate floor\n\n" + "top stars\n" + headline(summaries["top"], [t.name for t in real]))
if "total" in args.excess_scalings and main_scaling != "total":
    s_tot = pick(summaries["top"], scaling="total")
    s_tot = s_tot[s_tot["config"] == "2of3"].set_index("target")["r50_median"]
    head += "\n\n2 of 3 with the excess noise fixed in total-signal units ('total'): " + ", ".join(
        f"{t.name} {s_tot[t.name]:.0f} m" for t in real if t.name in s_tot.index)
theory = pick(summaries["top"], noise="theory")
theory = theory[theory["config"] == "2of3"]
head += "\n\n2 of 3, theory only [m]: " + ", ".join(f"{r.target} {r.r50_median:.0f}" for r in theory.itertuples())
head += "\n\ndummy stars, 2 of 3 (measured, theory) [m]:\n" + ", ".join(
    f"B={r.mag_B:g}: {r.r50_median:.0f} ({t.r50_median:.0f})" for r, t in zip(
        pick(summaries["dummy"]).query("config == '2of3'").sort_values("mag_B").itertuples(),
        pick(summaries["dummy"], noise="theory").query("config == '2of3'").sort_values("mag_B").itertuples()))
say("\n" + head)
with open(os.path.join(OUT, "headline.md"), "w") as f:
    f.write(head + "\n")
n_above = sum(int((t["r50_status"] == "above_grid").sum()) for t in tables.values())
n_below = sum(int((t["r50_status"] == "below_grid").sum()) for t in tables.values())
say(f"rows not reached within the grid: {n_above} above, {n_below} below")

# ───────────────────────────── search-bank check ────────────────────────────
# the 4-template bank of the real search (star of that size and temperature) against the perfect template at the star's R50
say("\nsearch-bank check ...")
instruments = {kind: Instrument(kind, filter_type=cfg.filter_type, site="LaPalma") for kind in ("MAGIC", "LST")}
first_run = [p for p in array_profiles if (p.night, p.observation) == (array_profiles[0].night, array_profiles[0].observation)]
rows = []
s2 = pick(summaries["top"]).query("config == '2of3'").set_index("target")["r50_median"]
for star in real:
    r50_star = float(s2.get(star.name, np.nan))
    if not np.isfinite(r50_star):
        continue
    cfg_b = replace(cfg, radii_m=np.array([r50_star]))
    bank = ms.bank_deltas(cfg_b, star=star)
    for p in first_run:
        kind = "LST" if p.telescope.startswith("LST") else "MAGIC"
        tpl = ms.template_family(star, cfg_b, instruments[kind])[float(r50_star)]
        psd = ms.noise_psd(p, star, cfg_b)
        gain = ms.correction_gain(cfg_b.correction, p.freq, psd, p.fs)
        perfect = ms.expected_snr(tpl, p.freq, psd, p.fs, gain)
        found = ms.bank_snr(tpl, bank, p.freq, psd, p.fs, gain)
        rows.append({"target": star.name, "mag_B": star.mag_B, "telescope": p.telescope, "radius_m": r50_star,
                     "snr_perfect": perfect, "snr_bank": found, "ratio": found / perfect})
bank_check = pd.DataFrame(rows)
bank_check.to_csv(os.path.join(OUT, "bank_check.csv"), index=False)
if len(bank_check):
    say(bank_check.groupby("target")["ratio"].agg(["min", "median", "max"]).round(3).to_string())

# ───────────────────────────── figures ──────────────────────────────────────
generic = None
if os.path.exists(args.generic):
    generic = pick(pd.read_csv(args.generic))
    generic = generic[generic["config"] == "2of3"].sort_values("mag_B")


def plot_r50(xcol, xlabel, fname, logx=False):
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    style(ax)
    d = pick(summaries["dummy"]).query("config == '2of3'").sort_values("mag_B")
    ax.plot(d[xcol], d["r50_median"], color=INK, lw=1.6, label=f"dummy stars, '{main_scaling}' excess")
    if len(args.excess_scalings) > 1:
        d2 = pick(summaries["dummy"], scaling=args.excess_scalings[-1]).query("config == '2of3'").sort_values("mag_B")
        ax.fill_between(d[xcol], d["r50_median"], d2["r50_median"], color=INK, alpha=0.12, lw=0,
                        label=f"up to '{args.excess_scalings[-1]}' excess")
    if "r50_lo_median" in d and d["r50_lo_median"].notna().any():
        ax.fill_between(d[xcol], d["r50_lo_median"], d["r50_hi_median"], color=INK, alpha=0.10, lw=0, label="calibration band")
    t = pick(summaries["dummy"], noise="theory").query("config == '2of3'").sort_values("mag_B")
    ax.plot(t[xcol], t["r50_median"], color=INK, ls="--", lw=1.2, label="dummy stars, theory only")
    if generic is not None and len(generic):
        g = generic.sort_values(xcol if xcol in generic else "mag_B")
        ax.plot(g[xcol], g["r50_median"], color=MUTED, ls=":", lw=1.4, label="generic 0.03 mas, 10 000 K stars")
    s = pick(summaries["top"]).query("config == '2of3'")
    ax.errorbar(s[xcol], s["r50_median"], yerr=[s["r50_median"] - s["r50_min"], s["r50_max"] - s["r50_median"]],
                ls="", marker="o", ms=6, color=ORANGE, ecolor=ORANGE, capsize=2.5, lw=1.2, zorder=5, label="top stars (run range)")
    for r in s.itertuples():
        ax.annotate(r.target, (getattr(r, xcol), r.r50_median), textcoords="offset points", xytext=(5, 5), fontsize=7, color=MUTED)
    if len(args.excess_scalings) > 1:
        s2_ = pick(summaries["top"], scaling=args.excess_scalings[-1]).query("config == '2of3'")
        ax.plot(s2_[xcol], s2_["r50_median"], ls="", marker="o", ms=6, markerfacecolor="none", color=ORANGE, zorder=5,
                label=f"top stars, '{args.excess_scalings[-1]}' excess")
    if logx:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel, fontsize=8, color=MUTED)
    ax.set_ylabel("R50, 2 of 3 [m]  (50 % detection)", fontsize=8, color=MUTED)
    if xcol == "star_mas":
        top_ax = ax.secondary_xaxis("top", functions=(lambda x: x * mas_to_rad * cfg.distance_au * AU_m,
                                                      lambda x: x / (mas_to_rad * cfg.distance_au * AU_m)))
        top_ax.set_xlabel(f"projected stellar radius at {cfg.distance_au:g} AU [m]", fontsize=8, color=MUTED)
        top_ax.tick_params(colors=MUTED, labelsize=8)
        for axis in (top_ax.xaxis, ax.xaxis):                  # plain numbers on the log axes
            axis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            axis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(f"Smallest detectable KBO radius, FAR {HEAD_FAR:g} / h, measured noise ({cfg.correction})", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=140)
    plt.close(fig)


plot_r50("mag_B", "B magnitude of the star", "r50_vs_mag.png")
plot_r50("star_mas", "angular radius of the star [mas]", "r50_vs_size.png", logx=True)

say(f"\nwrote {OUT}/")
