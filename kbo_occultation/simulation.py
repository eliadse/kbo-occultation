# kbo_occultation/simulation.py

import numpy as np

from .physics import (
    planck_photon,
    filter_transmission,
    fresnel_intensity_radial,
    AU_m,
    km_m,
    nm_m,
    mas_to_rad,
)


def compute_lightcurve(kbo, star, bandpass, grid, numerics):
    """
    Compute polychromatic occultation light curve.

    Returns
    -------
    x_km : array
    intensity : array
    """

    # ─── Derived quantities ─────────────────────────────

    D_m = kbo.distance_au * AU_m
    R_m = kbo.radius_km * km_m

    # Spatial grid
    x_km = np.linspace(-grid.x_max_km, grid.x_max_km, grid.n_x)
    x_m = x_km * km_m

    # Wavelength grid
    lambdas_nm = np.linspace(
        bandpass.lam_min_nm,
        bandpass.lam_max_nm,
        bandpass.n_lambda
    )
    lambdas_m = lambdas_nm * nm_m

    # Spectral weights
    spec_w = planck_photon(lambdas_m, star.temperature_K)
    filt_w = filter_transmission(
        lambdas_nm,
        bandpass.lam_min_nm,
        bandpass.lam_max_nm
    )

    weights = spec_w * filt_w
    weights /= weights.sum()

    # ─── Compute intensity ─────────────────────────────

    intensity_total = np.zeros_like(x_m)

    for lam, w in zip(lambdas_m, weights):
        I = fresnel_intensity_radial(
            x_m,
            R_m,
            D_m,
            lam,
            n_int=numerics.n_int
        )
        intensity_total += w * I

    return x_km, intensity_total
