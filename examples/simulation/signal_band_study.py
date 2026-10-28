"""
Signal-band study: where in frequency does the diffraction signal live?

1. Convergence: f-quantiles vs grid extent/resolution, wavelength sampling,
   stellar-disk sampling and radial-grid sampling. Acceptance: the two
   finest settings agree to < 2 % in f90 / f95.
2. Physical scan on the converged numerics: r*, R, b, Delta, v, star
   temperature, and the real MAGIC / LST instrument response.

Outputs (outputs/simulation/signal_band/): convergence.csv, scan.csv,
signal_band.png.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from kbo_occultation import BandpassConfig, GridConfig, NumericalConfig
from kbo_occultation.analysis.bands import signal_band
from kbo_occultation.simulation.instruments import Instrument

OUT_DIR = os.path.join("outputs", "simulation", "signal_band")
os.makedirs(OUT_DIR, exist_ok=True)

SITE = "La Palma"
V0 = 25e3
BAND0 = BandpassConfig(320, 450, 41)
GRID0 = GridConfig(30000, 12000)          # dx = 5 m
NUM0 = NumericalConfig(n_r_grid=6000, n_star_side=128)
BASE = dict(star_mas=0.03, radii_m=[100], distance_au=40.0, v_mps=V0, temperature_K=10000.0)
QUICK = "--quick" in sys.argv


NUM_SCAN = NumericalConfig(n_r_grid=6000, n_star_side=64)   # within ~1 % of 128 in f90, 4x faster
_numerics = NUM0


def run(**kw):
    args = dict(BASE, band=BAND0, grid=GRID0, numerics=_numerics)
    args.update(kw)
    return signal_band(**args)


# ─── 1. convergence ─────────────────────────────────────────────────
conv_rows = []
settings = [
    ("x_max", [dict(grid=GridConfig(10000, 4000)), dict(grid=GridConfig(30000, 12000)),
               dict(grid=GridConfig(60000, 24000))]),
    ("dx",    [dict(grid=GridConfig(30000, 6000)), dict(grid=GridConfig(30000, 12000)),
               dict(grid=GridConfig(30000, 24000))]),
    ("n_lambda", [dict(band=BandpassConfig(320, 450, n)) for n in (11, 41, 81)]),
    ("n_star_side", [dict(numerics=NumericalConfig(n_r_grid=6000, n_star_side=n)) for n in (16, 32, 64, 128, 256)]),
    ("n_r_grid", [dict(numerics=NumericalConfig(n_r_grid=n, n_star_side=128)) for n in (1500, 3000, 6000, 12000)]),
]
if QUICK:
    settings = settings[:3]
labels = {
    "x_max": ["10 km", "30 km", "60 km"], "dx": ["10 m", "5 m", "2.5 m"],
    "n_lambda": ["11", "41", "81"], "n_star_side": ["16", "32", "64", "128", "256"],
    "n_r_grid": ["1500", "3000", "6000", "12000"],
}
for name, variants in settings:
    for lab, kw in zip(labels[name], variants):
        for star in (0.01, 0.02, 0.03):
            r = run(star_mas=star, **kw)[0]
            conv_rows.append(dict(parameter=name, setting=lab, **{k: v for k, v in r.items()
                                  if k in ("star_mas", "radius_m") or k.endswith("_hz")}))
conv = pd.DataFrame(conv_rows)
conv.to_csv(os.path.join(OUT_DIR, "convergence.csv"), index=False)

print("\nConvergence: change between the two finest settings (%)")
ok = True
for name, group in conv.groupby("parameter", sort=False):
    settings_order = list(dict.fromkeys(group["setting"]))
    a, b = settings_order[-2], settings_order[-1]
    for star in (0.01, 0.02, 0.03):
        ra = group[(group.setting == a) & (group.star_mas == star)].iloc[0]
        rb = group[(group.setting == b) & (group.star_mas == star)].iloc[0]
        d90 = 100 * abs(rb.f90_hz / ra.f90_hz - 1)
        d95 = 100 * abs(rb.f95_hz / ra.f95_hz - 1)
        d99 = 100 * abs(rb.f99_hz / ra.f99_hz - 1)
        flag = "OK " if max(d90, d95) < 2 else "FAIL"
        ok &= flag == "OK "
        print(f"  {name:12s} {a:>6s}->{b:<6s} r*={star:.2f}  f90 {d90:5.2f}  f95 {d95:5.2f}  f99 {d99:5.2f}  {flag}")
print("ACCEPTANCE:", "PASS" if ok else "FAIL (see rows above)")

# ─── 2. physical scan ───────────────────────────────────────────────
_numerics = NUM_SCAN
scan = []


def add(tag, rows, **extra):
    for r in rows:
        scan.append(dict(scan=tag, **extra, **r))


STARS = np.round(np.arange(0.01, 0.0501, 0.005), 4)
RADII = [20, 50, 100, 200, 500, 1000]
for s in (STARS[::2] if QUICK else STARS):
    add("r_star_R", run(star_mas=float(s), radii_m=RADII))

REP_R = [50, 100, 300]
if not QUICK:
    for b in (0, 500, 1000, 2000):
        add("impact", run(radii_m=REP_R, impact_parameter_m=float(b)))
    for d in (30, 40, 50):
        add("distance", run(radii_m=REP_R, distance_au=float(d)))
    for T in (5000, 10000, 20000, 30000):
        add("temperature", run(radii_m=REP_R, temperature_K=float(T)))
    for tel in ("MAGIC", "LST"):
        inst = Instrument(telescope_type=tel, filter_type=None, site=SITE)
        add("response", run(radii_m=REP_R, response=inst.total_transmission), telescope=tel)
    # velocity: f scales exactly with v (f = f_x * v); evaluate directly
    base = run(radii_m=REP_R)
    for v in (5e3, 15e3, 25e3, 30e3):
        add("velocity", run(radii_m=REP_R, v_mps=v))

df = pd.DataFrame(scan)
df.to_csv(os.path.join(OUT_DIR, "scan.csv"), index=False)

# ─── figure ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
rs = df[df.scan == "r_star_R"]
for R, g in rs.groupby("radius_m"):
    axes[0].plot(g.star_mas, g.f90_hz, "-o", ms=3, label=f"R={R:g} m")
    axes[1].plot(g.star_mas, g.f_low_hz, "-o", ms=3, label=f"R={R:g} m")
axes[0].set(xlabel="stellar radius (mas)", ylabel="f90 (Hz)", title="Upper band edge (90 % of power)")
axes[1].set(xlabel="stellar radius (mas)", ylabel="f05 (Hz)", title="Lower band edge (5 % of power)")
axes[1].set_yscale("log")
for ax in axes[:2]:
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
g = rs[rs.radius_m == 100]
for q in ("f50_hz", "f90_hz", "f95_hz", "f99_hz"):
    axes[2].plot(g.star_mas, g[q], "-o", ms=3, label=q[:-3])
axes[2].set(xlabel="stellar radius (mas)", ylabel="f (Hz)", title="Quantiles, R=100 m (f99 is tail-dominated)")
axes[2].legend(fontsize=8); axes[2].grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "signal_band.png"), dpi=180)

pd.set_option("display.width", 200)
print("\nBand (Hz) at v = 25 km/s, Delta = 40 AU, b = 0, R = 100 m, MAGIC-like 320-450 nm, T = 10 000 K")
print(g[["star_mas", "f_low_hz", "f50_hz", "f90_hz", "f95_hz", "f99_hz", "bandwidth_hz", "chord_time_s"]]
      .round(2).to_string(index=False))
