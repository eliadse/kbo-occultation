"""
``kbo-night``: take an observation night from raw files to an SNR table.

    kbo-night list                       nights in the log and what is ready / missing
    kbo-night prepare 2026-10-20         write/complete the night in the log (SIMBAD, CC, AMC,
                                         channels), build .npz caches + conditions, list what's missing
    kbo-night conditions 2026-10-20      observing conditions per observation
    kbo-night snr 2026-10-20 --plot      measured vs theoretical SNR -> outputs/snr/

NIGHT arguments accept several nights; none (or "all") means every night.
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
    _show(pd.DataFrame(status).set_index("night").join(table[["telescopes", "n_stars", "dark_run"]]))
    print(f"\nlog: {log.path}\n`kbo-night prepare NIGHT` lists the missing items of a night.")


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
    _show(conditions_table(log, args.nights), args.csv)


def cmd_snr(args):
    import matplotlib
    matplotlib.use("Agg")
    from .analysis.snr import plot_snr_vs_magnitude, snr_report

    log = _load(args)
    nights = [n.name for n in log.select(args.nights)]
    results = snr_report(log, nights, n_bins=args.bins, detrend_window_s=args.detrend_window)
    tag = "all" if len(nights) == len(log.nights) else "_".join(nights)
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


def build_parser():
    p = argparse.ArgumentParser(prog="kbo-night", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--log", default=DEFAULT_LOG_PATH, help="observation_log.yaml to use")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("list", help="nights in the log and their status")
    s.set_defaults(func=cmd_list)

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
    s.set_defaults(func=cmd_conditions)

    s = sub.add_parser("snr", help="measured vs theoretical SNR")
    s.add_argument("nights", nargs="*")
    s.add_argument("--bins", type=int, nargs="+", default=[10, 200], help="samples per bin")
    s.add_argument("--detrend-window", type=float, default=10.0, help="baseline window for detr [s]")
    s.add_argument("--out", default="outputs/snr", help="output directory")
    s.add_argument("--plot", action="store_true", help="also plot SNR vs magnitude per bin size")
    s.set_defaults(func=cmd_snr)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    try:
        args.func(args)
    except (ObservationLogError, FileNotFoundError) as e:
        sys.exit(f"kbo-night: {e}")


if __name__ == "__main__":
    main()
