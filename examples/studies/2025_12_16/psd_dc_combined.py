# examples/studies/2025_12_16/dc_combine.py
"""
Apply the fast/slow-DC frequency-domain combination (dc_combine.py) to every
star of one night in observation_log.yaml and check whether it actually
removes the correlated sub-few-Hz noise excess found in
examples/studies/2025_12_16/inspect_real_noise.py.

PSD of the fractional-flux fluctuations, raw vs. DC-combined, one row per
star and one column per channel -- if the correlated component really is
coming from whatever the slow DC already sees, the combined PSD should show
a reduced sub-few-Hz excess.

Usage:
    python examples/studies/2025_12_16/dc_combine.py [NIGHT] [--channels A B C] [--cut-hz 1.0]
"""

import argparse
import os

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import welch

from kbo_occultation.config import DC_CUT_HZ
from kbo_occultation.observations.log import load_observation_log
from kbo_occultation.observations.dc_combine import combine_lightcurve_with_dc

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("night", nargs="?", default="2025-12-16", help="night in observation_log.yaml (YYYY-MM-DD)")
parser.add_argument("--channels", nargs="+", default=None, help="channels to use (default: the night's channels)")
parser.add_argument("--cut-hz", type=float, default=None,
                    help="fast/DC crossover frequency (default: per-telescope DC_CUT_HZ from analysis.snr)")
args = parser.parse_args()

log = load_observation_log()
night = log.select([args.night])[0]
channels = args.channels or night.channels

if night.dc_report is None or not os.path.exists(night.dc_report):
    raise SystemExit(f"night {night.name}: DC report not found ({night.dc_report})")

OUT_DIR = os.path.join("outputs", "studies", night.name)
os.makedirs(OUT_DIR, exist_ok=True)

stars = []
for obs in night.stars.values():
    if not obs.has_cache:
        print(f"skip {obs.name}: no .npz cache (run `kbo-night prepare {night.name}`)")
        continue
    stars.append(obs)
stars.sort(key=lambda o: o.mag_B if o.mag_B is not None else 99)
if not stars:
    raise SystemExit(f"night {night.name}: no star observations with a .npz cache")


def psd(lc):
    norm = lc.signal / np.mean(lc.signal) - 1.0
    nperseg = min(2**16, len(norm))
    return welch(norm, fs=1.0 / lc.dt, window="hann", nperseg=nperseg, noverlap=nperseg // 2, scaling="density")


# ─── PSD: raw vs. combined ──────────────────────────────────────────────
fig, axes = plt.subplots(len(stars), len(channels), figsize=(9 * len(channels), 2.6 * len(stars)),
                         sharex=True, squeeze=False)
for i, obs in enumerate(stars):
    for j, ch in enumerate(channels):
        ax = axes[i, j]
        tel = night.telescope(ch)
        cut_hz = args.cut_hz if args.cut_hz is not None else DC_CUT_HZ[tel.dc_key]
        raw = obs.lightcurve(ch, good_time=True)
        curves = [("raw", raw)]
        try:
            curves.append(("DC-combined", combine_lightcurve_with_dc(
                raw, night.dc_report, frequency_cut_hz=cut_hz, telescope=tel.dc_key, pixel=0)))
        except (ValueError, KeyError) as err:
            print(f"{obs.name} / {tel.name}: no DC combination ({err})")
        for label, lc in curves:
            f, Pxx = psd(lc)
            ax.plot(f[1:], Pxx[1:], lw=1.1, label=label)
        ax.axvline(cut_hz, color="k", ls="--", lw=1.0, alpha=0.6, label="cut")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(r"PSD of $\delta I/I$")
        ax.set_title(f"{obs.name} (B={obs.mag_B}, V={obs.mag_V}) -- {tel.name} (ch {ch})")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3, which="both")
for ax in axes[-1]:
    ax.set_xlabel("Frequency (Hz)")
fig.tight_layout()
out_name = f"noise_psd_dc_combined_{night.name}.png"
fig.savefig(os.path.join(OUT_DIR, out_name), dpi=200)
plt.show()

print(f"\nSaved: {os.path.join(OUT_DIR, out_name)}")
