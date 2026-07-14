"""
Diffraction pattern, power spectral density, and cumulative power for KBOs
of different radii, star-light only (no instrument response applied).

Three stacked panels, one column per KBO radius worth of curves overlaid:

1. Diffraction pattern I(x), noiseless, b=0, on the native spatial grid.
2. Power spectral density vs spatial frequency (m^-1), on the same grid
   (spatial_power_spectrum operates directly on x_m, before any
   time conversion, so units stay spatial throughout).
3. Cumulative power vs spatial frequency (m^-1) -- fraction of total
   spectral power contained below each frequency.
"""

import numpy as np
import matplotlib.pyplot as plt

from kbo_occultation import GridConfig, NumericalConfig, BandpassConfig, run_parameter_sweep
from kbo_occultation.physics import AU_m, nm_m
from kbo_occultation.detectability import (
    spatial_to_time, resample_to_cadence, sigma_from_instrument,
    spatial_power_spectrum, nyquist_frequency,
)

STAR_TEMPERATURE_K = 30000
SHADOW_VELOCITY_MPS = 25e3
STAR_ANGULAR_RADIUS_MAS = 0.03
DISTANCE_AU = 40.0
RADII_M = [50, 100, 200, 500]
THRESHOLD = 0.99

BAND = BandpassConfig(lam_min_nm=400, lam_max_nm=500, n_lambda=41)

grid = GridConfig(x_max_m=6000, n_x=100)
numerics = NumericalConfig(n_r_grid=2000)

# Fresnel scale F = sqrt(lambda * D / 2), evaluated at the passband's mean
# wavelength -- used only to put a second, achromatic-ish set of axes
# (and the KBO radii) in Fresnel units alongside the physical m / m^-1 ones.
LAM_MEAN_M = 0.5 * (BAND.lam_min_nm + BAND.lam_max_nm) * nm_m
D_M = DISTANCE_AU * AU_m
FRESNEL_SCALE_M = np.sqrt(LAM_MEAN_M * D_M / 2.0)

# responses=None -> each band falls back to its own trapezoidal
# filter_transmission weighting, i.e. just the star's spectrum over the
# passband, with no instrument response folded in.
points = run_parameter_sweep(
    radii_m=RADII_M,
    distances_au=[DISTANCE_AU],
    impact_parameters_m=[0.0],
    star_angular_radii_mas=[STAR_ANGULAR_RADIUS_MAS],
    star_temperature_K=STAR_TEMPERATURE_K,
    bands=[BAND],
    grid=grid,
    numerics=numerics,
    keep_curves=True,
)

#curves = {p.radius_m: (p.x_m, p.intensity) for p in points}
# On Bickerton 2009 they plot I/<I>
curves = {p.radius_m: (p.x_m, p.intensity) for p in points}

fig, (ax_curve, ax_psd, ax_cum) = plt.subplots(3, 1, figsize=(8, 10), sharex=False)

for radius_m in sorted(curves):
    x_m, I = curves[radius_m]
    freq_per_m, power = spatial_power_spectrum(x_m, I)
    cumulative = np.cumsum(power) / np.sum(power)
    rho = radius_m / FRESNEL_SCALE_M
    # Getting the Nyquist frequency

    f_nyq_hz, dt_max_s = nyquist_frequency(x_m, I, SHADOW_VELOCITY_MPS)
    print("F) Nyquist freq., max dt: ", f_nyq_hz, dt_max_s * 1e3, "ms")

    label = f"R={radius_m} m ($\\rho$={rho:.2f} F)"
    ax_curve.plot(x_m, I, label=label, lw=1.3)
    ax_psd.plot(freq_per_m, power, label=label, lw=1.3)
    ax_cum.plot(freq_per_m, cumulative, label=label, lw=1.3)

ax_curve.set_xlabel("x (m)")
ax_curve.set_ylabel("Relative intensity")
ax_curve.set_title(f"Diffraction pattern, b=0, star size={STAR_ANGULAR_RADIUS_MAS} mas, "
                    f"{BAND.lam_min_nm:g}-{BAND.lam_max_nm:g} nm")
ax_curve.legend(fontsize=8)
ax_curve.grid(alpha=0.3)
ax_curve_top = ax_curve.secondary_xaxis("top", functions=(lambda x: x / FRESNEL_SCALE_M,
                                                           lambda r: r * FRESNEL_SCALE_M))
ax_curve_top.set_xlabel("Position in shadow (Fsu)")

#ax_psd.set_yscale("log")
ax_psd.set_xlim(0, 0.0015)
ax_psd.set_xlabel("Wave Number (m$^{-1}$)")
ax_psd.set_ylabel("Normalised PSD")
ax_psd.legend(fontsize=8)
ax_psd.grid(alpha=0.3)
ax_psd_top = ax_psd.secondary_xaxis("top", functions=(lambda f: f * FRESNEL_SCALE_M,
                                                       lambda fr: fr / FRESNEL_SCALE_M))
ax_psd_top.set_xlabel("Wave Number ($\mathrm{Fsu^{-1}}$)")

ax_cum.set_xlim(0, 0.0015)
ax_cum.set_xlabel("Wave number (m$^{-1}$)")
ax_cum.set_ylabel("Cum power")
ax_cum.set_xlim(ax_psd.get_xlim())
ax_cum.legend(fontsize=8)
ax_cum.grid(alpha=0.3)
ax_cum.axhline(THRESHOLD, color="black", ls="--", lw=1.0, alpha=0.7)
ax_cum_top = ax_cum.secondary_xaxis("top", functions=(lambda f: f * FRESNEL_SCALE_M,
                                                       lambda fr: fr / FRESNEL_SCALE_M))
ax_cum_top.set_xlabel("Wave Number ($\mathrm{Fsu^{-1}}$)")

fig.tight_layout()
#fig.savefig("diffraction_psd_cumulative_BKO09_check.png", dpi=200)
fig.savefig(f"diffraction_psd_cumulative_{STAR_ANGULAR_RADIUS_MAS}.png", dpi=200)
plt.show()
