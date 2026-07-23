# examples/array_mc_benchmark.py
"""
Full-array injection/recovery Monte Carlo on a benchmark-star
observation, and the resulting upper limit on the cumulative KBO
surface density.

Flow:
1. Load all three channels (C=MAGIC-1, A=MAGIC-2, B=LST-1) from one
   stat-binary file.
2. Inject the same simulated occultation into every channel (with the
   ms-scale inter-telescope offsets) at random times, for a grid of KBO
   radii, and recover each injection per telescope + in coincidence
   (>= MIN_NTEL telescopes within the coincidence tolerance).
3. Turn the coincidence recovery fraction into an efficiency curve
   eps(r), and eps(r) + the exposure (live time, shadow velocity,
   distance) into a 95% CL upper limit N(>r0) [deg^-2].

Note on the 2025-12-16 data used here by default: only MAGIC-2
(channel A) had a real star signal that night, so the A/B/C channels
are not three independent views of the same star -- treat the output as
a pipeline exercise. With the 2026 three-telescope benchmark runs this
script applies as-is.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import LightCurve
from kbo_occultation.config import (BandpassConfig, GridConfig, NumericalConfig,
                                    ORM_ARRAY, StarConfig)
from kbo_occultation.kinematics import ShadowKinematics, shadow_velocity_opposition
from kbo_occultation.injection_batch import (efficiency_curve,
                                             run_array_injection_monte_carlo,
                                             to_dataframe)
from kbo_occultation.filtering import highpass_lightcurve
from kbo_occultation.upper_limits import SurveyExposure, cumulative_density_upper_limit

DATA_DIR = f"{PACKAGE_DATA}/observations"
DATA_FILE = f"{DATA_DIR}/Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_HD-17316_10002_20002_20251216T223319.npz"
STAR_NAME = "HD-17316"

DISTANCE_AU = 43.0
STAR = StarConfig(temperature_K=30000, angular_radius_mas=0.03)
BAND = BandpassConfig(300.0, 650.0, 25)
GRID = GridConfig(x_max_m=4000, n_x=800)
NUMERICS = NumericalConfig(n_r_grid=3000)

RADII_M = [100, 200, 300, 500]
N_TRIALS_PER_RADIUS = 20
PAD_S = 20.0
MIN_NTEL = 2
HIGHPASS_CUTOFF_HZ = 3.0
SEED = 2026

# Injection/exposure velocity: the opposition value the benchmark run
# emulates, so eps(r) carries over to on-ecliptic data (see
# kinematics module docstring for the split convention).
KIN = ShadowKinematics(v_rel_mps=shadow_velocity_opposition(DISTANCE_AU))

method_specs = {
    "raw": lambda lc: lc,
    "highpass": lambda lc: highpass_lightcurve(lc, HIGHPASS_CUTOFF_HZ),
}

lcs = LightCurve.from_preprocessed_all(DATA_FILE)

trials, profiles = run_array_injection_monte_carlo(
    lcs, ORM_ARRAY,
    radii_m=RADII_M, n_trials_per_radius=N_TRIALS_PER_RADIUS,
    grid=GRID, numerics=NUMERICS, star=STAR, band=BAND, response=None,
    kinematics=KIN, distance_au=DISTANCE_AU,
    method_specs=method_specs, pad_s=PAD_S, min_ntel=MIN_NTEL, seed=SEED,
    # highpass reshapes the event's own fine structure (an LTI filter),
    # so the search/veto must compare against the *filtered* template --
    # without this, large filtered events fail the chi2 shape veto
    # against the unfiltered template and eps(r) collapses.
    reshape_template_methods=("highpass",),
    checkpoint_path="array_mc_benchmark_trials.csv",
)
to_dataframe(profiles).to_csv("array_mc_benchmark_noise_profiles.csv", index=False)

live_time_s = float(lcs["A"].time[-1] - lcs["A"].time[0])
radii_grid = np.linspace(min(RADII_M), max(RADII_M), 100)

fig, (ax_eff, ax_ul) = plt.subplots(1, 2, figsize=(12, 5))
for method in method_specs:
    eff = efficiency_curve(trials, method, min_ntel=MIN_NTEL)
    exposure = SurveyExposure(
        star_name=STAR_NAME, live_time_s=live_time_s,
        shadow_velocity_mps=KIN.v_rel_mps, distance_au=DISTANCE_AU,
        efficiency=eff,
    )
    ul = cumulative_density_upper_limit([exposure], radii_grid)

    ax_eff.plot(eff.radii_m, eff.efficiency, "o-", label=method)
    ax_ul.plot(ul.radii_m, ul.n_upper_deg2, label=method)
    print(f"[{method}] eps(r): " +
          ", ".join(f"{r:.0f}m={e:.2f}" for r, e in zip(eff.radii_m, eff.efficiency)))

ax_eff.set_xlabel("KBO radius [m]")
ax_eff.set_ylabel(f"coincidence efficiency (>= {MIN_NTEL} telescopes)")
ax_eff.set_ylim(0, 1.05)
ax_eff.legend()
ax_ul.set_xlabel("r0 [m]")
ax_ul.set_ylabel(r"N(>r0) 95% CL upper limit [deg$^{-2}$]")
ax_ul.set_yscale("log")
ax_ul.legend()
fig.suptitle(f"{STAR_NAME}, {live_time_s:.0f}s live time, single observation")
fig.tight_layout()
fig.savefig("array_mc_benchmark.png", dpi=150)
print("wrote array_mc_benchmark.png, array_mc_benchmark_trials.csv, "
      "array_mc_benchmark_noise_profiles.csv")
