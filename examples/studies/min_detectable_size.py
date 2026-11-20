# examples/studies/min_detectable_size.py
"""
Smallest KBO that MAGIC-1 + MAGIC-2 + LST-1 can detect (kbo_occultation.analysis.min_size, steps 1-3).

R50 is the radius at which the detection probability is 50 % with >= 2 of the 3 telescopes in coincidence at a
coincidence false-alarm rate of 1 per hour (0.1 and 0.01 per hour alongside). The noise is the measured one
of the logged star runs (the red noise that the plain photon budget misses), carried to other magnitudes; "theory"
is the ideal system of the noise model. The expected SNRs carry the calibration of
``min_size_calibration.py`` (``SNR_CALIBRATION``) and its run-to-run spread is reported as a band.

Targets: a generic grid B = 7 ... 9 (angular radius 0.03 mas * 10^(-0.2 (B - 7)), 39 AU, opposition, b = 0),
and the stars of ``--catalogue``. The LST-1 infinity-focus scenario multiplies the star light of LST-1 by 2; how the
excess noise follows ("star": fixed in star units, "total": fixed in total-signal units) is not known from the data,
so both are run.

Besides R50 the detection width W(R) = 2 * integral of the detection probability over the impact parameter b is
computed (``min_size.detection_width``) and turned into the observing time per detected event with the size
distribution of Bickerton et al. 2009 (``rates.observing_time_per_size``, ``rates.detected_event_rate``).

Usage:
    python examples/studies/min_detectable_size.py [NIGHT ...] [--min-channels 3] [--fars 1 0.1 0.01]
        [--mags 7 9 0.25] [--radii 10 1000 41] [--focus 1 2] [--excess-scalings star total]
        [--catalogue stars.csv --columns name=Name,mag_B=Bmag,angular_radius_mas=theta] [--no-b-average]
        [--check-numerics] [--out outputs/min_size]
Smoke test:
    python examples/studies/min_detectable_size.py --mags 7 9 1 --radii 10 1000 15 --fars 1 --no-b-average \
        --out outputs/min_size/smoke
"""
import argparse
import os
import sys
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u

from kbo_occultation.analysis import min_size as ms
from kbo_occultation.analysis import rates
from kbo_occultation.observations.log import load_observation_log

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("nights", nargs="*", default=["2025-12-16", "2026-03-08", "2026-08-19"],
                help="nights of the noise profiles and rate curves (the array statistics use runs with --min-channels)")
ap.add_argument("--min-channels", type=int, default=3, help="the headline uses the runs with this many channels")
ap.add_argument("--correction", default="whitened", help="raw | whitened | butter_<fc>")
ap.add_argument("--fars", nargs="+", type=float, default=[1.0, 0.1, 0.01], help="coincidence FARs [1/h], the first is the headline")
ap.add_argument("--mags", nargs=3, type=float, default=[7.0, 9.0, 0.25], metavar=("MIN", "MAX", "STEP"))
ap.add_argument("--radii", nargs=3, type=float, default=[10.0, 1000.0, 41], metavar=("MIN", "MAX", "N"),
                help="KBO radius grid [m], geometric")
ap.add_argument("--focus", nargs="+", type=float, default=[1.0, 2.0], help="LST-1 star-light gains (2 = infinity focus)")
ap.add_argument("--excess-scalings", nargs="+", default=["star", "total"])
ap.add_argument("--width-mags", nargs="+", type=float, default=[7.0, 8.0, 9.0], help="magnitudes of the b-averaged widths")
ap.add_argument("--no-b-average", action="store_true", help="skip the impact-parameter averaging (slow step)")
ap.add_argument("--d-break-km", type=float, default=80.0)
ap.add_argument("--q-s", type=float, default=4.0, help="small-object slope of the size distribution")
ap.add_argument("--catalogue", default=None, help="CSV of target stars")
ap.add_argument("--columns", default=None, help="catalogue column map: name=Name,mag_B=Bmag,angular_radius_mas=...")
ap.add_argument("--check-numerics", action="store_true", help="n_star_side / n_r_grid convergence and search-bank mismatch")
ap.add_argument("--stationary", action="store_true", help="one PSD per run instead of one per time block")
ap.add_argument("--profiles", default=None, help="noise profiles .npz (default outputs/min_size/profiles_<nights>.npz)")
ap.add_argument("--candidates", default=os.path.join("outputs", "noise_removal", "real_candidates.csv"))
ap.add_argument("--calibration", default=os.path.join("outputs", "min_size", "calibration"),
                help="directory of min_size_calibration.py (snr_calibration_per_run.csv gives the band)")
ap.add_argument("--out", default=os.path.join("outputs", "min_size"))
ap.add_argument("--include-excluded", action="store_true")
args = ap.parse_args()

OUT = args.out
os.makedirs(OUT, exist_ok=True)
COLORS = {"MAGIC-2": "#2a78d6", "MAGIC-1": "#eb6834", "LST-1": "#1baf7a"}     # validated categorical slots 1-3
MARKERS = {"MAGIC-2": "o", "MAGIC-1": "s", "LST-1": "^", "2of3": "D", "3of3": "v"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
ARRAY_COLORS = {"2of3": INK, "3of3": MUTED}
D_BREAK = args.d_break_km * u.km


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


def line_style(config):
    color = COLORS.get(config, ARRAY_COLORS.get(config, INK))
    return dict(color=color, marker=MARKERS.get(config, "o"), ms=4.5, lw=1.5,
                ls="--" if config == "3of3" else "-")


log = load_observation_log()
radii = np.geomspace(args.radii[0], args.radii[1], int(args.radii[2]))
cfg = ms.MinSizeConfig(radii_m=radii, correction=args.correction, coinc_far_per_hr=args.fars[0],
                       extra_fars=tuple(args.fars[1:]), time_resolved=not args.stationary).calibrated()
say(f"correction {cfg.correction}, snr_scale {cfg.snr_scale:.3f}, scatter {cfg.snr_scatter:.3f}, "
    f"{'time resolved' if cfg.time_resolved else 'stationary'} noise, FAR {args.fars}")
HEAD_FAR = cfg.coinc_far_per_hr
tol = cfg.tolerance_s()

# ───────────────────────────── noise, thresholds, targets ───────────────────
prof_path = args.profiles or os.path.join("outputs", "min_size", "profiles_" + "_".join(args.nights) + ".npz")
profiles = ms.load_or_measure_profiles(log, args.nights, cfg, prof_path, include_excluded=args.include_excluded,
                                       progress=say)
runs = ms.calibration_runs(log, args.nights, profiles, include_excluded=args.include_excluded,
                           min_channels=args.min_channels)
if not runs:
    sys.exit(f"no run with >= {args.min_channels} channels in {args.nights}")
array_profiles = [p for _, _, group in runs for p in group]
groups = {(n.name, o.name): g for n, o, g in runs}
say(f"{len(runs)} runs: " + ", ".join(f"{n.name}/{o.name}" for n, o, _ in runs))

curves = ms.rate_curves(args.candidates, log, ms.csv_method_name(cfg.correction), args.nights,
                        include_excluded=args.include_excluded)
thresholds = ms.thresholds_table(curves, tol, args.fars)
thresholds.to_csv(os.path.join(OUT, "thresholds.csv"), index=False)
say("thresholds (per-telescope SNR):")
say(thresholds.round(2).to_string(index=False))

band = None
per_run_csv = os.path.join(args.calibration, "snr_calibration_per_run.csv")
if os.path.exists(per_run_csv):
    try:
        band = ms.scale_quantiles(per_run_csv, cfg.correction)
        say(f"calibration band (16-84 % of the per-run scale): {band[0]:.2f} - {band[1]:.2f} (used: {cfg.snr_scale:.2f})")
    except ValueError as e:
        say(f"no calibration band: {e}")

targets = ms.generic_targets(np.arange(args.mags[0], args.mags[1] + 1e-9, args.mags[2]))

# ───────────────────────────── R50 ──────────────────────────────────────────
say("R50 table ...")
table, eff = ms.min_size_table(targets, array_profiles, cfg, curves, focus_gains=args.focus,
                               excess_scalings=args.excess_scalings, scale_band=band, return_efficiency=True,
                               progress=say)
table.to_csv(os.path.join(OUT, "r50_table.csv"), index=False)
summary = ms.summarize(table)
summary.to_csv(os.path.join(OUT, "r50_summary.csv"), index=False)
head = eff[(eff["far_per_hr"] == HEAD_FAR) & (eff["lst_focus"] == 1.0)
           & eff["excess_scaling"].isin([cfg.excess_scaling, "none"])]
head.groupby(["target", "mag_B", "noise", "config", "radius_m"], as_index=False)["eff"].median() \
    .to_csv(os.path.join(OUT, "efficiency.csv"), index=False)

main_scaling = args.excess_scalings[0]


def pick(df, noise="measured", focus=1.0, far=None, scaling=None):
    far = HEAD_FAR if far is None else far
    q = df[(df["noise"] == noise) & (df["lst_focus"] == focus) & (df["far_per_hr"] == far)]
    if noise == "measured":
        q = q[q["excess_scaling"] == (scaling or main_scaling)]
    return q


def headline(config_order):
    s = pick(summary)
    s = s[s["config"].isin(config_order)]
    lines = ["| configuration | " + " | ".join(f"B = {m:g}" for m in (7.0, 8.0, 9.0)) + " |", "|---|---|---|---|"]
    for cfg_name in config_order:
        cells = []
        for m in (7.0, 8.0, 9.0):
            r = s[(s["config"] == cfg_name) & np.isclose(s["mag_B"], m)]
            if r.empty:
                cells.append("-")
                continue
            r = r.iloc[0]
            lo, hi = r.get("r50_lo_median", np.nan), r.get("r50_hi_median", np.nan)
            text = f"{r.r50_median:.0f} ({r.r50_min:.0f}-{r.r50_max:.0f})"
            if np.isfinite(lo) and np.isfinite(hi):
                text += f" [{lo:.0f}-{hi:.0f}]"
            cells.append(text + ("*" if r.thr_at_floor else ""))
        lines.append(f"| {cfg_name} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


tel_names = sorted({p.telescope for p in array_profiles})
configs_present = tel_names + [c for c in ("2of3", "3of3") if c in set(table["config"])]
say(f"\nR50 [m], measured noise, excess '{main_scaling}', FAR {HEAD_FAR:g}/h: median (run range) [calibration band]; "
    "* = threshold at the candidate floor")
say(headline(configs_present))
theory = summary[(summary["noise"] == "theory") & (summary["lst_focus"] == 1.0) & (summary["far_per_hr"] == HEAD_FAR)
                 & (summary["config"] == "2of3") & summary["mag_B"].isin([7.0, 8.0, 9.0])]
say("theory-only 2 of 3 [m]: " + ", ".join(f"B={r.mag_B:g}: {r.r50_median:.0f}" for r in theory.itertuples()))
unreached = int((table["r50_status"] == "above_grid").sum()), int((table["r50_status"] == "below_grid").sum())
say(f"rows not reached within the grid: {unreached[0]} above, {unreached[1]} below")

# ───────────────────────────── catalogue ────────────────────────────────────
if args.catalogue:
    stars = ms.load_target_catalogue(args.catalogue, ms.parse_column_map(args.columns))
    say(f"catalogue: {len(stars)} stars")
    cat = ms.min_size_table(stars, array_profiles, replace(cfg, extra_fars=()), curves, noise_models=("measured",),
                            focus_gains=args.focus, excess_scalings=args.excess_scalings, scale_band=band, progress=say)
    cat_summary = ms.summarize(cat)
    cat_summary.to_csv(os.path.join(OUT, "catalogue_r50.csv"), index=False)
    say(cat_summary[(cat_summary["config"] == "2of3") & (cat_summary["lst_focus"] == 1.0)]
        [["target", "mag_B", "excess_scaling", "r50_median", "r50_min", "r50_max"]].round(0).to_string(index=False))

# ───────────────────────────── detection widths and observing time ──────────
width_summary = pd.DataFrame()
if not args.no_b_average:
    from kbo_occultation.simulation.instruments import Instrument

    say("\nimpact-parameter averaged detection widths ...")
    radii_w = radii[::2] if len(radii) > 21 else radii
    cfg_w = replace(cfg, radii_m=radii_w)
    instruments = {kind: Instrument(kind, filter_type=cfg.filter_type, site="LaPalma") for kind in ("MAGIC", "LST")}
    rows = []
    for mag in args.width_mags:
        star = ms.generic_targets([mag])[0]
        templates = {kind: ms.template_family(star, cfg_w, ins) for kind, ins in instruments.items()}
        signals = {kind: ms.width_signals(star, cfg_w, ins) for kind, ins in instruments.items()}
        say(f"  B = {mag:g}: signals done")
        for (night, obs), group in groups.items():
            configs = [(p.telescope, [p], 1) for p in group]
            if len(group) >= 2:
                configs.append((f"2of{len(group)}", group, 2))
            for noise in ("measured", "theory"):
                for label, sel, k in configs:
                    if noise == "theory" and not label.startswith("2of"):
                        continue
                    chans = tuple(p.channel for p in sel)
                    thr, _ = ms.array_threshold(curves, chans, k, HEAD_FAR, tol)
                    w = ms.detection_width(star, sel, cfg_w, templates, signals, k=k, threshold=thr,
                                           theory_only=(noise == "theory"))
                    rows.append(w.assign(mag_B=mag, night=night, observation=obs, noise=noise, config=label,
                                         threshold=thr))
    widths = pd.concat(rows, ignore_index=True)
    widths.to_csv(os.path.join(OUT, "detection_width.csv"), index=False)
    width_summary = (widths.groupby(["noise", "config", "mag_B", "radius_m"], as_index=False)
                     .agg(h_m=("h_m", "median"), width_m=("width_m", "median"), eff_eff=("eff_eff", "median"),
                          p_b0=("p_b0", "median"), p_tail=("p_tail", "max"), n_runs=("width_m", "size")))
    width_summary["label"] = width_summary["config"] + " " + width_summary["noise"]
    width_summary.to_csv(os.path.join(OUT, "detection_width_summary.csv"), index=False)
    if (width_summary["p_tail"] > 0.02).any():
        bad = width_summary[width_summary["p_tail"] > 0.02]
        say(f"warning: {len(bad)} cells have P > 0.02 at the outermost impact parameter (W is cut short)")

    velocity = cfg.velocity_mps() * u.m / u.s
    distance = (cfg.distance_au * u.au).to(u.m)
    obs_rows = []
    for (label, mag), g in width_summary.groupby(["label", "mag_B"]):
        g = g.sort_values("radius_m")
        t_size = rates.observing_time_per_size(g["radius_m"].to_numpy(), g["width_m"].to_numpy(), velocity, distance,
                                               D_BREAK, args.q_s).to(u.hr).value
        mu = rates.detected_event_rate(g["radius_m"].to_numpy(), g["width_m"].to_numpy(), velocity, distance,
                                       D_BREAK, args.q_s)
        for r, t in zip(g["radius_m"], t_size):
            obs_rows.append({"label": label, "mag_B": mag, "radius_m": r, "hours_per_event": t,
                             "hours_per_event_integrated": float((1.0 / mu).to(u.hr).value)})
    obs_time = pd.DataFrame(obs_rows)
    obs_time.to_csv(os.path.join(OUT, "observing_time.csv"), index=False)
    integrated = obs_time.drop_duplicates(["label", "mag_B"]).pivot(index="label", columns="mag_B",
                                                                    values="hours_per_event_integrated")
    say(f"\nobserving time per detected event, whole size distribution (q_s = {args.q_s:g}, break {args.d_break_km:g} km) [h]:")
    say(integrated.round(1).to_string())

    # the old shot-noise model for the array of three MAGIC (tests/observing_time_vs_kbo_size.py)
    shot = []
    for mag in args.width_mags:
        star = ms.generic_targets([mag])[0]
        ins = instruments["MAGIC"]
        thr = 1.0 / float(ins.signal_to_noise_ratio(mag, 100 * u.ms, ntels=3))
        for r, (x_m, intensity) in ms.injection_patterns(star, cfg_w, ins, radii_w).items():
            dev = np.flatnonzero(np.abs(intensity - 1.0) > thr)
            dev = dev[(dev >= 15) & (dev < len(x_m) - 15)]
            shot.append({"mag_B": mag, "radius_m": r, "width_m": float(x_m[dev[-1]] - x_m[dev[0]]) if dev.size else 0.0})
    shot = pd.DataFrame(shot)
    shot["hours_per_event"] = [float(rates.observing_time_per_size([r], [w], velocity, distance, D_BREAK, args.q_s)
                                     .to(u.hr).value[0]) for r, w in zip(shot["radius_m"], shot["width_m"])]
    shot.to_csv(os.path.join(OUT, "observing_time_shot_noise_model.csv"), index=False)

# ───────────────────────────── numerics ─────────────────────────────────────
if args.check_numerics:
    say("\nnumerics checks ...")
    from kbo_occultation.simulation.instruments import Instrument

    cfg_n = replace(cfg, radii_m=np.geomspace(max(radii[0], 20.0), radii[-1], 8))
    group0 = next(iter(groups.values()))
    instruments = {kind: Instrument(kind, filter_type=cfg.filter_type, site="LaPalma") for kind in ("MAGIC", "LST")}
    bank = ms.bank_deltas(cfg_n)
    rows = []
    for mag in (7.0, 9.0):
        star = ms.generic_targets([mag])[0]
        base = {kind: ms.template_family(star, cfg_n, ins) for kind, ins in instruments.items()}
        variants = {"n_star_side 64 -> 128": replace(cfg_n, n_star_side=128),
                    "n_r_grid 4000 -> 8000": replace(cfg_n, n_r_grid=8000)}
        for name, cfg_v in variants.items():
            alt = {kind: ms.template_family(star, cfg_v, ins) for kind, ins in instruments.items()}
            for p in group0:
                kind = "LST" if p.telescope.startswith("LST") else "MAGIC"
                a = ms.per_telescope_snr(star, p, base[kind], cfg_n)
                b = ms.per_telescope_snr(star, p, alt[kind], cfg_n)
                a, b = (np.nanmedian(v, axis=0) if v.ndim == 2 else v for v in (a, b))
                rows += [{"check": name, "mag_B": mag, "telescope": p.telescope, "radius_m": r, "a": x, "b": y,
                          "rel_diff": y / x - 1.0} for r, x, y in zip(cfg_n.radii_m, a, b)]
        bank_star = ms.bank_deltas(cfg_n, star=star)         # the 4-template bank of this star's own size
        for p in group0:                                       # the 4-template banks against a perfect template
            kind = "LST" if p.telescope.startswith("LST") else "MAGIC"
            f_noise = ms.noise_psd(p, star, cfg_n)
            gain = ms.correction_gain(cfg_n.correction, p.freq, f_noise, p.fs)
            for r in cfg_n.radii_m:
                perfect = ms.expected_snr(base[kind][float(r)], p.freq, f_noise, p.fs, gain)
                for name, bk in (("per-star bank vs perfect template", bank_star),
                                 ("fixed 0.03 mas bank vs perfect template", bank)):
                    found = ms.bank_snr(base[kind][float(r)], bk, p.freq, f_noise, p.fs, gain)
                    rows.append({"check": name, "mag_B": mag, "telescope": p.telescope, "radius_m": r,
                                 "a": perfect, "b": found, "rel_diff": found / perfect - 1.0})
    num = pd.DataFrame(rows)
    num.to_csv(os.path.join(OUT, "numerics.csv"), index=False)
    say(num.groupby("check")["rel_diff"].agg(median="median", min="min", max="max",
                                              max_abs=lambda x: float(np.abs(x).max())).round(4).to_string())

# ───────────────────────────── figures ──────────────────────────────────────
mags_plot = sorted(set(summary["mag_B"]))
fars = args.fars
fig, axes = plt.subplots(1, len(fars), figsize=(4.6 * len(fars), 4.2), sharey=True, squeeze=False)
for ax, far in zip(axes[0], fars):
    style(ax)
    s = pick(summary, far=far)
    for cfg_name in configs_present:
        g = s[s["config"] == cfg_name].sort_values("mag_B")
        if g.empty:
            continue
        kw = line_style(cfg_name)
        ax.plot(g["mag_B"], g["r50_median"], label=cfg_name if cfg_name not in ("2of3", "3of3") else
                {"2of3": "2 of 3", "3of3": "3 of 3"}[cfg_name], **kw)
        floor = g[g["thr_at_floor"].astype(bool)]
        ax.plot(floor["mag_B"], floor["r50_median"], ls="", marker=kw["marker"], ms=kw["ms"] + 1.5,
                markerfacecolor="#fcfcfb", markeredgecolor=kw["color"], zorder=4)
        if cfg_name == "2of3":
            ax.fill_between(g["mag_B"], g["r50_min"], g["r50_max"], color=INK, alpha=0.10, lw=0, label="run range")
            if "r50_lo_median" in g and g["r50_lo_median"].notna().any():
                ax.fill_between(g["mag_B"], g["r50_lo_median"], g["r50_hi_median"], color=INK, alpha=0.18, lw=0,
                                label="calibration band")
    t = pick(summary, noise="theory", far=far)
    t = t[t["config"] == "2of3"].sort_values("mag_B")
    if len(t):
        ax.plot(t["mag_B"], t["r50_median"], color=INK, ls=":", lw=1.5, label="2 of 3, theory only")
    ax.set_xlabel("B magnitude of the star", fontsize=8, color=MUTED)
    ax.set_title(f"coincidence FAR {far:g} / h", fontsize=9, color=INK)
axes[0, 0].set_ylabel("R50 [m]  (50 % detection)", fontsize=8, color=MUTED)
handles, labels = axes[0, 0].get_legend_handles_labels()
for ax in axes[0, 1:]:                                   # series that only some panels have
    for h, lab in zip(*ax.get_legend_handles_labels()):
        if lab not in labels:
            handles.append(h)
            labels.append(lab)
fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 5), fontsize=7.5, frameon=False)
fig.suptitle(f"Smallest detectable KBO radius, measured noise ({cfg.correction}, excess '{main_scaling}')\n"
             "open marker: threshold at the candidate floor (the FAR is then only an upper bound)",
             fontsize=9, color=INK)
fig.tight_layout(rect=[0, 0.1, 1, 0.93])
fig.savefig(os.path.join(OUT, "r50_vs_mag.png"), dpi=140)
plt.close(fig)

if len(args.focus) > 1 and any(p.telescope.startswith("LST") for p in array_profiles):
    cfgs = [c for c in ("LST-1", "2of3") if c in configs_present]
    fig, axes = plt.subplots(1, len(cfgs), figsize=(4.8 * len(cfgs), 4.2), sharey=True, squeeze=False)
    for ax, cfg_name in zip(axes[0], cfgs):
        style(ax)
        color = line_style(cfg_name)["color"]
        for focus, ls in ((1.0, "-"), (2.0, "--")):
            lines = {sc: pick(summary, focus=focus, scaling=sc).query("config == @cfg_name").sort_values("mag_B")
                     for sc in args.excess_scalings}
            first = lines[args.excess_scalings[0]]
            ax.plot(first["mag_B"], first["r50_median"], color=color, ls=ls, lw=1.6,
                    label=f"focus x{focus:g}, excess '{args.excess_scalings[0]}'")
            if len(lines) > 1:
                last = lines[args.excess_scalings[-1]]
                ax.fill_between(first["mag_B"], first["r50_median"], last["r50_median"], color=color,
                                alpha=0.14 if focus == 1.0 else 0.30, lw=0,
                                label=f"focus x{focus:g}: up to '{args.excess_scalings[-1]}'")
        ax.set_xlabel("B magnitude of the star", fontsize=8, color=MUTED)
        ax.set_title("2 of 3" if cfg_name == "2of3" else cfg_name, fontsize=9, color=INK)
        ax.legend(fontsize=7, frameon=False)
    axes[0, 0].set_ylabel("R50 [m]", fontsize=8, color=MUTED)
    fig.suptitle(f"LST-1 at infinity focus (star light x2), FAR {HEAD_FAR:g} / h: how the excess noise follows is not known",
                 fontsize=9, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(OUT, "r50_lst_focus.png"), dpi=140)
    plt.close(fig)

if len(width_summary):
    cols = sorted(set(width_summary["mag_B"]))
    fig, axes = plt.subplots(1, len(cols), figsize=(4.4 * len(cols), 4.0), sharey=True, squeeze=False)
    for ax, mag in zip(axes[0], cols):
        style(ax)
        for noise, ls in (("measured", "-"), ("theory", ":")):
            g = width_summary[(width_summary["mag_B"] == mag) & (width_summary["noise"] == noise)
                              & width_summary["config"].str.startswith("2of")].sort_values("radius_m")
            if g.empty:
                continue
            ax.plot(g["radius_m"], g["p_b0"], color=INK, ls=ls, lw=1.3, label=f"b = 0, {noise}")
            if noise == "measured":
                ax.plot(g["radius_m"], g["eff_eff"], color=COLORS["LST-1"], marker="o", ms=3.5, lw=1.5,
                        label="W / H(R), averaged over b")
        ax.axhline(1.0, color=MUTED, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("KBO radius [m]", fontsize=8, color=MUTED)
        ax.set_title(f"B = {mag:g}", fontsize=9, color=INK)
    axes[0, 0].set_ylabel("efficiency (2 of 3)", fontsize=8, color=MUTED)
    axes[0, -1].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "efficiency_curves.png"), dpi=140)
    plt.close(fig)

    shot = pd.read_csv(os.path.join(OUT, "observing_time_shot_noise_model.csv"))
    fig, axes = plt.subplots(1, len(cols), figsize=(4.4 * len(cols), 4.0), sharey=True, squeeze=False)
    for ax, mag in zip(axes[0], cols):
        style(ax)
        for label, kw in ((f"{g} measured", line_style(g)) for g in [*tel_names, "2of3"] if g in configs_present):
            g = obs_time[(obs_time["label"] == label) & (obs_time["mag_B"] == mag)].sort_values("radius_m")
            if len(g):
                ax.plot(2.0 * g["radius_m"], g["hours_per_event"], label=label.replace("measured", "").strip()
                        .replace("2of3", "2 of 3"), **kw)
        g = obs_time[(obs_time["label"] == "2of3 theory") & (obs_time["mag_B"] == mag)].sort_values("radius_m")
        if len(g):
            ax.plot(2.0 * g["radius_m"], g["hours_per_event"], color=INK, ls=":", lw=1.5, label="2 of 3, theory only")
        g = shot[shot["mag_B"] == mag].sort_values("radius_m")
        g = g[np.isfinite(g["hours_per_event"])]
        ax.plot(2.0 * g["radius_m"], g["hours_per_event"], color=MUTED, ls="--", lw=1.2, label="shot-noise model, 3 MAGIC")
        for hours, text in ((24, "1 day"), (168, "1 week"), (730, "1 month"), (8766, "1 year")):
            ax.axhline(hours, color=GRID, lw=0.8)
            ax.text(2.0 * radii[0] * 1.05, hours, text, fontsize=7, color=MUTED, va="bottom", ha="left")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(3.0, 1e8)                       # below a few hours / above 10^4 years nothing is informative
        ax.set_xlabel("KBO diameter [m]", fontsize=8, color=MUTED)
        ax.set_title(f"B = {mag:g}", fontsize=9, color=INK)
    axes[0, 0].set_ylabel("observing time per detected event [h]", fontsize=8, color=MUTED)
    axes[0, -1].legend(fontsize=7, frameon=False, loc="upper right")
    fig.suptitle(f"Observing time for one detected event of objects of that size (q_s = {args.q_s:g}, "
                 f"{cfg.distance_au:g} AU, opposition, FAR {HEAD_FAR:g} / h)", fontsize=9, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(OUT, "observing_time.png"), dpi=140)
    plt.close(fig)

say(f"\nwrote {OUT}/")
