# examples/studies/noise_removal_comparison.py
"""
Which noise-removal method keeps the most detection power?

Efficiency eps(R, mag) at a fixed real-data false-alarm rate (see
kbo_occultation.analysis.noise_removal): Butterworth cutoff scan, PSD-whitened
matched filter, DC detrend window scan, fixed DC FFT swap, DC despike as a veto.

Usage:
    python examples/studies/noise_removal_comparison.py [NIGHT ...] [--far 1.0] [--trials 20]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kbo_occultation.analysis.noise_removal import ComparisonConfig, run_comparison, summarize
from kbo_occultation.observations.log import load_observation_log

OUT_DIR = os.path.join("outputs", "noise_removal")

ap = argparse.ArgumentParser()
ap.add_argument("nights", nargs="*")
ap.add_argument("--far", type=float, default=1.0, help="real-data false alarms per telescope-hour")
ap.add_argument("--trials", type=int, default=20, help="injections per run, channel, radius, method")
ap.add_argument("--methods", nargs="+", default=None, help="only these methods (default: all)")
ap.add_argument("--out", default=OUT_DIR, help="output directory")
ap.add_argument("--include-excluded", action="store_true")
args = ap.parse_args()

os.makedirs(OUT_DIR, exist_ok=True)
OUT_DIR = args.out
os.makedirs(OUT_DIR, exist_ok=True)
cfg = ComparisonConfig(far_per_hr=args.far, n_trials=args.trials, methods=args.methods)
thr, eff, cand = run_comparison(load_observation_log(), args.nights or None, cfg,
                                include_excluded=args.include_excluded, progress=print)
thr.to_csv(f"{OUT_DIR}/thresholds.csv", index=False)
eff.to_csv(f"{OUT_DIR}/efficiency.csv", index=False)
cand.to_csv(f"{OUT_DIR}/real_candidates.csv", index=False)
table = summarize(eff)
table.to_csv(f"{OUT_DIR}/efficiency_summary.csv")
print(thr.to_string(index=False))
print(table.round(2).to_string())

fig, (ax_r, ax_m) = plt.subplots(1, 2, figsize=(13, 5))
for method, g in eff.groupby("method"):
    by_r = g.groupby("radius_m")["efficiency"].mean()
    ax_r.plot(by_r.index, by_r.values, "o-", label=method)
    if g["mag_B"].notna().any():
        by_m = g.groupby("mag_B")["efficiency"].mean()
        ax_m.plot(by_m.index, by_m.values, "o-", label=method)
ax_r.set_xscale("log"); ax_r.set_xlabel("KBO radius [m]"); ax_r.set_ylabel(f"efficiency at FAR = {args.far:g} / telescope-hour")
ax_m.set_xlabel("star B magnitude"); ax_m.set_ylabel("efficiency (mean over radii)")
for ax in (ax_r, ax_m):
    ax.set_ylim(-0.02, 1.05); ax.grid(alpha=0.3)
ax_r.legend(fontsize=7)
fig.tight_layout()
fig.savefig(f"{OUT_DIR}/efficiency.png", dpi=150)
print(f"wrote {OUT_DIR}/")
