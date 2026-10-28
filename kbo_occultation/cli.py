"""
``kbo-night``: take an observation night from raw files to an SNR table.

    kbo-night list                       nights in the log and what is ready / missing
    kbo-night missing 2026-10-20         what is still missing for the night (read-only)
    kbo-night prepare 2026-10-20         write/complete the night in the log (SIMBAD, CC, AMC,
                                         channels), build .npz caches + conditions, list what's missing
    kbo-night conditions 2026-10-20      observing conditions per observation
    kbo-night snr 2026-10-20 --plot      measured vs theoretical SNR -> outputs/snr/
    kbo-night mimic 2025-12-16           candidates, Gaussian surrogates and false-alarm rates -> outputs/mimic/

NIGHT arguments accept several nights; none (or "all") means every night.
Nights / observations with ``exclude:`` in the log are skipped by ``conditions``
and ``snr`` unless named (nights) or ``--include-excluded`` is given.
``--log`` selects another observation_log.yaml.
"""

import argparse
import os
import sys
import warnings

import pandas as pd

from .observations.log import DEFAULT_LOG_PATH, ObservationLogError, load_observation_log, nights_table, validate


def _load(args):
    return load_observation_log(args.log)


def _show(df: pd.DataFrame, csv=None):
    with pd.option_context("display.width", 250, "display.max_columns", 60, "display.precision", 3):
        print(df.to_string())
    if csv:
        os.makedirs(os.path.dirname(os.path.abspath(csv)), exist_ok=True)
        df.to_csv(csv)
        print(f"\nwritten {csv}")


def _select(log, args):
    """log.select for the command's nights, reporting on stderr what was excluded."""
    if not args.include_excluded:
        for item in log.excluded(args.nights):
            print(f"excluded {item} (--include-excluded to use it)", file=sys.stderr)
    return log.select(args.nights, args.include_excluded)


def cmd_list(args):
    log = _load(args)
    table = nights_table(log)
    status = []
    for night in log.nights.values():
        obs = list(night.observations.values())
        n_cached = sum(o.has_cache for o in obs)
        has_cond = bool(night.conditions_dir and os.path.isdir(night.conditions_dir))
        status.append({"night": night.name, "runs": len(obs), "cached": n_cached,
                       "conditions": "yes" if has_cond else ("no" if night.cc_reports else "n/a"),
                       "missing": len(validate(night))})
    _show(pd.DataFrame(status).set_index("night").join(table[["telescopes", "n_stars", "dark_run", "exclude"]]))
    print(f"\nlog: {log.path}\n`kbo-night missing NIGHT` lists the missing items of a night (read-only);"
          f"\n`kbo-night prepare NIGHT` fills in what it can and lists what's still missing.")


def cmd_missing(args):
    for night in _load(args).select(args.nights, include_excluded=True):
        problems = validate(night)
        print(f"── {night.name}: missing")
        for p in problems:
            print(f"  {p}")
        if not problems:
            print("  nothing")


def cmd_prepare(args):
    from .observations.prepare import prepare_night

    dates = args.nights
    if not dates or dates == ["all"]:
        dates = list(_load(args).nights)
    for date in dates:
        prepare_night(date, log_path=args.log, force=args.force, refresh=args.refresh,
                      dry_run=args.dry_run, pad_s=args.pad)


def cmd_conditions(args):
    from .observations.cc_reports import conditions_table

    log = _load(args)
    nights = [n.name for n in _select(log, args)]
    _show(conditions_table(log, nights, args.include_excluded), args.csv)


def cmd_snr(args):
    import matplotlib
    matplotlib.use("Agg")
    from .analysis.snr import plot_snr_vs_magnitude, snr_report

    log = _load(args)
    nights = [n.name for n in _select(log, args)]
    results = snr_report(log, nights, n_bins=args.bins, detrend_window_s=args.detrend_window,
                         include_excluded=args.include_excluded)
    tag = "all" if not args.nights or args.nights == ["all"] else "_".join(nights)
    os.makedirs(args.out, exist_ok=True)
    csv = os.path.join(args.out, f"snr_{tag}.csv")
    results.to_csv(csv)
    cols = [c for c in ("detr", "mean_fast", "theo_obs_B", "ratio_detr_B", "airmass", "mirror_throughput")
            if c in results]
    _show(results[cols])
    print(f"\nwritten {csv} (all columns; see `kbo_occultation.analysis.snr` for their meaning)")
    if args.plot:
        for n_bin in args.bins:
            fig = plot_snr_vs_magnitude(results, log, n_bin)
            png = os.path.join(args.out, f"snr_{tag}_{n_bin}bin.png")
            fig.savefig(png, dpi=120)
            print(f"written {png}")


def cmd_mimic(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .analysis.mimic import (MimicConfig, far_table, flag_enrichment_table, run_mimic, tail_table)

    log = _load(args)
    nights = [n.name for n in _select(log, args)]
    cfg = MimicConfig(corrections=tuple(args.corrections), snr_threshold=args.snr_threshold,
                      n_surrogates=args.surrogates, n_shifts=args.shifts, bin_samples=args.bin_samples,
                      star_mas=args.star_mas, distance_au=args.distance_au, seed=args.seed,
                      sigma_window_s=args.sigma_window)
    res = run_mimic(log, nights, cfg, include_excluded=args.include_excluded,
                    progress=lambda m: print(m, file=sys.stderr))
    tag = "all" if not args.nights or args.nights == ["all"] else "_".join(nights)
    out = os.path.join(args.out, tag)
    os.makedirs(out, exist_ok=True)
    far, tails, enrich = far_table(res), tail_table(res), flag_enrichment_table(res)
    for name, df in (("candidates", res.candidates), ("events", res.events), ("live", res.live),
                     ("far", far), ("tails", tails), ("flag_enrichment", enrich),
                     ("calibration", res.calibration)):
        df.to_csv(os.path.join(out, f"{name}.csv"), index=False)
    if len(far):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        for ax, scope in zip(axes, ("single", "coincident")):
            for (corr, kind), g in far[far.scope == scope].groupby(["correction", "kind"]):
                ax.errorbar(g.snr_min, g.rate_per_hr.clip(lower=1e-3), g.rate_err, marker="o", ms=3,
                            ls={"real": "-", "surrogate": "--", "timeslide": ":"}[kind],
                            label=f"{corr} / {kind}")
            ax.set(yscale="log", xlabel="SNR >=", ylabel="events / hour", title=scope)
            ax.grid(alpha=0.3)
        axes[0].legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "far.png"), dpi=130)
    print("\nFalse-alarm rate (events / hour) at SNR >= 5")
    if len(far):
        _show(far[far.snr_min == 5.0].set_index(["correction", "scope", "kind"])[["count", "hours", "rate_per_hr", "rate_err"]])
    print(f"\nwritten {out}/ (candidates, events, live, far, tails, flag_enrichment, calibration .csv, far.png)")


def build_parser():
    p = argparse.ArgumentParser(prog="kbo-night", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--log", default=DEFAULT_LOG_PATH, help="observation_log.yaml to use")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("list", help="nights in the log and their status")
    s.set_defaults(func=cmd_list)

    s = sub.add_parser("missing", help="what is still missing for nights in the log (read-only)")
    s.add_argument("nights", nargs="*")
    s.set_defaults(func=cmd_missing)

    s = sub.add_parser("prepare", help="create/complete nights in the log, build caches and conditions")
    s.add_argument("nights", nargs="*", help="night dates YYYY-MM-DD (new or already in the log); none = all")
    s.add_argument("--force", action="store_true", help="rebuild existing caches / conditions")
    s.add_argument("--refresh", action="store_true",
                   help="overwrite automatic fields already set in the log (manual fields are never touched)")
    s.add_argument("--dry-run", action="store_true", help="print the log changes without writing them")
    s.add_argument("--pad", type=float, default=1800.0, help="CC window padding around the runs [s]")
    s.set_defaults(func=cmd_prepare)

    s = sub.add_parser("conditions", help="observing conditions per observation")
    s.add_argument("nights", nargs="*")
    s.add_argument("--csv", help="also write the table to this CSV")
    s.add_argument("--include-excluded", action="store_true", help="also use nights/observations marked exclude")
    s.set_defaults(func=cmd_conditions)

    s = sub.add_parser("snr", help="measured vs theoretical SNR")
    s.add_argument("nights", nargs="*")
    s.add_argument("--bins", type=int, nargs="+", default=[10, 200], help="samples per bin")
    s.add_argument("--detrend-window", type=float, default=10.0, help="baseline window for detr [s]")
    s.add_argument("--out", default="outputs/snr", help="output directory")
    s.add_argument("--plot", action="store_true", help="also plot SNR vs magnitude per bin size")
    s.add_argument("--include-excluded", action="store_true", help="also use nights/observations marked exclude")
    s.set_defaults(func=cmd_snr)

    s = sub.add_parser("mimic", help="how often noise looks like an occultation (candidates, surrogates, FAR)")
    s.add_argument("nights", nargs="*")
    s.add_argument("--night", dest="night_opt", action="append", help="night date (same as the positional)")
    s.add_argument("--corrections", nargs="+", default=["raw", "dc_detrend", "dc_despike"])
    s.add_argument("--snr-threshold", type=float, default=3.0)
    s.add_argument("--surrogates", type=int, default=20, help="Gaussian surrogates per run and correction")
    s.add_argument("--shifts", type=int, default=200, help="time-slide realizations per run")
    s.add_argument("--bin-samples", type=int, default=10, help="samples averaged before the search")
    s.add_argument("--star-mas", type=float, default=0.03, help="stellar angular radius of the templates")
    s.add_argument("--distance-au", type=float, default=40.0)
    s.add_argument("--sigma-window", type=float, default=None,
                   help="window [s] of the empirical filter-output sigma (default: detrend window)")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--out", default="outputs/mimic", help="output directory")
    s.add_argument("--include-excluded", action="store_true", help="also use nights/observations marked exclude")
    s.set_defaults(func=cmd_mimic)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if getattr(args, "night_opt", None):
        args.nights = list(args.nights) + args.night_opt
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    try:
        args.func(args)
    except (ObservationLogError, FileNotFoundError) as e:
        sys.exit(f"kbo-night: {e}")


if __name__ == "__main__":
    main()
