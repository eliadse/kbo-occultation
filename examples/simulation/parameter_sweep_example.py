"""
Example: sweep KBO radius, impact parameter, and stellar size across two
filters, then flag which combinations would be detectable at a given
(placeholder, flat) noise level.

Swap `flat_sigma` below for your actual per-observation noise estimate
(photon + read + scintillation + systematic, combined in quadrature) once
that's wired up -- detectability.matched_filter_snr just needs a sigma in
fractional-flux units, it doesn't care where it came from.
"""

import matplotlib.pyplot as plt
import numpy as np

from kbo_occultation import (
    BandpassConfig,
    GridConfig,
    NumericalConfig,
    run_parameter_sweep,
    to_records,
    matched_filter_snr,
    is_detectable,
)

grid = GridConfig(x_max_m=3000, n_x=400)
numerics = NumericalConfig(n_r_grid=2500)

points = run_parameter_sweep(
    radii_m=[50, 100, 200, 500],
    distances_au=[40.0],
    impact_parameters_m=[0.0, 20.0, 50.0, 100.0, 500.0],
    star_angular_radii_mas=[0.01, 0.02, 0.03, 0.04],
    star_temperature_K=15000,
    bands=[BandpassConfig(400, 430, 10), BandpassConfig(550, 700, 10)],
    grid=grid,
    numerics=numerics,
    keep_curves=True
)

# --- flat placeholder noise level; replace with your real noise model ---
flat_sigma = 0.03

for p in points:
    p_snr = matched_filter_snr(p.intensity, sigma=flat_sigma)
    p.snr = p_snr  # stash it on the object for convenience below
    p.detectable = p_snr >= 5.0

records = to_records(points)
print(f"{len(records)} sweep points computed")
print(f"detectable @ SNR>=5, sigma={flat_sigma}: "
      f"{sum(p.detectable for p in points)} / {len(points)}")

# --- quick detectability map: radius vs impact parameter, one filter,
#     smallest star size ---
band_label = "400-430nm"
star_mas = 0.01
subset = [
    p for p in points
    if p.band_label == band_label and p.star_angular_radius_mas == star_mas
]
radii = sorted(set(p.radius_m for p in subset))
impacts = sorted(set(p.impact_parameter_m for p in subset))
snr_grid = np.zeros((len(radii), len(impacts)))
for p in subset:
    i = radii.index(p.radius_m)
    j = impacts.index(p.impact_parameter_m)
    snr_grid[i, j] = p.snr

fig, ax = plt.subplots()
im = ax.imshow(
    snr_grid, origin="lower", aspect="auto",
    extent=[min(impacts), max(impacts), min(radii), max(radii)],
)
ax.set_xlabel("Impact parameter (m)")
ax.set_ylabel("KBO radius (m)")
ax.set_title(f"Matched-filter SNR, {band_label}, star={star_mas} mas, sigma={flat_sigma}")
fig.colorbar(im, label="SNR")
plt.show()
