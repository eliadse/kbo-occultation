"""
Visual comparison of the diffraction light curves themselves (not just
summary SNR) for MAGIC vs LST, filtered vs unfiltered. Companion to
magic_lst_filter_comparison.py -- keep the parameters below in sync with
that script if you change one.

Two figures:

1. lightcurves_noiseless.png -- noiseless light curves for every KBO
   radius in RADII_M, b=0, one panel per instrument+filter config. This
   is where you can actually *see* what the SNR numbers are summarizing:
   broadband (unfiltered) responses integrate over a wide range of
   Fresnel fringe spacings (which scale with wavelength), so the fringes
   from different wavelengths fall out of phase with each other and wash
   out into one smooth, wide dip. The narrowband SII response keeps a
   much narrower wavelength range, so the fringes stay in phase and the
   diffraction structure (secondary fringes on either side of the main
   dip) survives much more clearly.

2. lightcurves_with_noise.png -- one representative radius, resampled to
   the real sampling cadence, with the noiseless curve and a single
   Gaussian noise realization (at each config's actual sigma from
   sigma_from_instrument) overlaid -- so you can sanity-check the SNR
   numbers against what the event would actually look like in the data,
   not just trust the summary statistic.
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from kbo_occultation import GridConfig, NumericalConfig, BandpassConfig, run_parameter_sweep
from kbo_occultation.config import standard_sampling, standard_sample_duration
from kbo_occultation.instruments import Instrument
from kbo_occultation.detectability import spatial_to_time, resample_to_cadence, sigma_from_instrument

# ─── Keep in sync with magic_lst_filter_comparison.py ───────────────
STAR_TEMPERATURE_K = 30000
STAR_MAGNITUDE = 8.0
SHADOW_VELOCITY_MPS = 25e3
SITE = "La Palma"
SAMPLING_DT_S = 20 * standard_sampling * standard_sample_duration * 1e-9

CONFIGS = [
    ("MAGIC", None),
    ("MAGIC", "MAGIC_SII"),
    ("LST", None),
    ("LST", "MAGIC_SII"),
]

RADII_M = [50, 100]
STAR_ANGULAR_RADII_MAS = [0.03]
DISTANCE_AU = 40.0
REPRESENTATIVE_RADIUS_M = 80  # used for the noisy single-event figure; must be in RADII_M

grid = GridConfig(x_max_m=4000, n_x=800)
numerics = NumericalConfig(n_r_grid=3000)

# label -> {radius_m: (t_s, intensity)}, all at b=0, full spatial resolution
curves = {}
sigmas = {}

for tel, filt in CONFIGS:
    inst = Instrument(telescope_type=tel, filter_type=filt, site=SITE)
    l_min, l_max, n = inst.response_band()
    band = BandpassConfig(l_min, l_max, n)
    response = inst.total_transmission_filter if filt else inst.total_transmission
    label = f"{tel}{'+SII' if filt else ' (unfiltered)'}"
    band_label = f"{band.lam_min_nm:g}-{band.lam_max_nm:g}nm"

    points = run_parameter_sweep(
        radii_m=RADII_M,
        distances_au=[DISTANCE_AU],
        impact_parameters_m=[0.0],
        star_angular_radii_mas=STAR_ANGULAR_RADII_MAS,
        star_temperature_K=STAR_TEMPERATURE_K,
        bands=[band],
        grid=grid,
        numerics=numerics,
        responses={band_label: response},
        keep_curves=True,
    )

    sigmas[label] = sigma_from_instrument(inst, STAR_MAGNITUDE, SAMPLING_DT_S * u.s)
    curves[label] = {
        p.radius_m: (spatial_to_time(p.x_m, SHADOW_VELOCITY_MPS), p.intensity)
        for p in points
    }

    print(f"{label:16s} band={band.lam_min_nm:5.1f}-{band.lam_max_nm:5.1f}nm  sigma={sigmas[label]:.4f}")

# ─── Figure 1: noiseless curves, all radii, one panel per config ────
fig1, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
for ax, (label, per_radius) in zip(axes.flat, curves.items()):
    for radius_m, (t_s, I) in sorted(per_radius.items()):
        ax.plot(t_s * 1e3, I, label=f"R={radius_m} m", lw=1.2)
    ax.set_title(label)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Relative intensity")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
fig1.suptitle("Noiseless diffraction light curves, b=0 (full spatial-grid resolution)")
fig1.tight_layout()
fig1.savefig("lightcurves_noiseless.png", dpi=200)

# ─── Figure 1: noiseless curves, all radii, one panel per config ────
fig1, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
for ax, (label, per_radius) in zip(axes.flat, curves.items()):
    for radius_m, (t_s, I) in sorted(per_radius.items()):
        # TODO add the fourier transform to plot and then get the Nyquist frequency
        # TODO I think I should get the frequency before doing the "spatial_to_time"
        ax.plot(t_s * 1e3, I, label=f"R={radius_m} m", lw=1.2)
    ax.set_title(label)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Relative intensity")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
fig1.suptitle("Noiseless diffraction light curves, b=0 (full spatial-grid resolution)")
fig1.tight_layout()
fig1.savefig("lightcurves_noiseless.png", dpi=200)


# ─── Figure 2: one radius, noiseless + a noise realization, all configs ──
rng = np.random.default_rng(0)
fig2, ax2 = plt.subplots(figsize=(8, 5))
offset = 0.0
for label, per_radius in curves.items():
    t_s, I = per_radius[REPRESENTATIVE_RADIUS_M]
    t_binned, I_binned = resample_to_cadence(t_s, I, SAMPLING_DT_S)
    sigma = sigmas[label]
    noisy = I_binned + rng.normal(0, sigma, size=I_binned.shape)

    ax2.plot(t_binned * 1e3, noisy + offset, ".", ms=4, alpha=0.6, color=f"C{list(curves).index(label)}")
    ax2.plot(t_binned * 1e3, I_binned + offset, "-", lw=1.5, color=f"C{list(curves).index(label)}",
              label=f"{label} (sigma={sigma:.3f})")
    offset += 0.4  # stack traces vertically so they don't overlap

ax2.set_xlabel("Time (ms)")
ax2.set_ylabel("Relative intensity (traces offset for clarity)")
ax2.set_title(f"R={REPRESENTATIVE_RADIUS_M} m, b=0 -- noiseless curve + one noise realization, real cadence")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)
fig2.tight_layout()
fig2.savefig("lightcurves_with_noise.png", dpi=200)

plt.show()
