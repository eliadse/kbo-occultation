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

def apply_stellar_disk(x_m, intensity, star_radius_m, n_star_side):
    """
    Convolve intensity with a uniform stellar disk.
    """

    # Build 2D grid of stellar offsets
    R = star_radius_m

    offsets = np.linspace(-R, R, n_star_side)
    dx, dy = np.meshgrid(offsets, offsets)

    mask = dx**2 + dy**2 <= R**2

    dx = dx[mask]
    dy = dy[mask]

    convolved = np.zeros_like(intensity)

    for shift in dx:
        shifted = np.interp(
            x_m + shift,
            x_m,
            intensity,
            left=1.0,
            right=1.0
        )
        convolved += shifted

    convolved /= len(dx)

    return convolved

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
    
    #if bandpass.lam_min_nm != bandpass.lam_max_nm:
    #    weights /= weights.sum()
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

    # projected stellar radius
    r_star_m = star.angular_radius_mas * mas_to_rad * D_m

    if numerics.n_star_side > 1:
        intensity_total = apply_stellar_disk(
            x_m,
            intensity_total,
            r_star_m,
            numerics.n_star_side
        )
    return x_km, intensity_total

#def simulate_poly_point(x_m, b_m, R_m, D_m, lambdas_m, weights, N_int=800):
def simulate_poly_point(kbo, star, bandpass, grid, numerics):
    
    # KBO parameters
    D_m = kbo.distance_au * AU_m
    R_m = kbo.radius_km * km_m
    b_m = kbo.impact_parameter_km * km_m
    
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

    N_int = numerics.n_int

    """ Monochromatic and polychromatic point source"""
    r_obs = np.sqrt(x_m**2 + b_m**2)
    total = np.zeros(len(x_m))
    for lam_m, w in zip(lambdas_m, weights):
        if w < 1e-12:
            continue
        total += w * fresnel_intensity_radial(r_obs, R_m, D_m, lam_m, N_int)
    return x_km, total