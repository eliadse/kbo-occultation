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
from .instruments import *

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

def apply_stellar_disk_2d(x_m, intensity_radial, r_grid_m, star_radius_m, impact_parameter_m, n_star_side):
    """
    Exact 2D convolution with a uniform stellar disk.

    Parameters
    ----------
    x_m : 1D array
        Chord positions
    intensity_radial : 1D array
        I(r) evaluated on r_grid_m
    r_grid_m : 1D array
        Radial grid corresponding to intensity_radial
    star_radius_m : float
        Projected stellar radius
    impact_parameter_m: float
        Distance of the KBO plane from the star center
    n_star_side : int
        Resolution of stellar disk grid

    Returns
    -------
    convolved_intensity : 1D array
    """

    # --- Build stellar disk sampling ---
    offsets = np.linspace(-star_radius_m, star_radius_m, n_star_side)
    dx, dy = np.meshgrid(offsets, offsets)

    mask = dx**2 + dy**2 <= star_radius_m**2

    dx = dx[mask]
    dy = dy[mask]

    # --- Interpolator for radial intensity ---
    def interp_I(r):
        return np.interp(r, r_grid_m, intensity_radial, left=1.0, right=1.0)

    # --- Convolution ---
    convolved = np.zeros_like(x_m)

    for i, x in enumerate(x_m):
        r = np.sqrt((x + dx)**2 + (impact_parameter_m + dy)**2)
        convolved[i] = interp_I(r).mean()

    return convolved


def compute_lightcurve(kbo, star, bandpass, grid, numerics, SII=False):
    """
    Compute polychromatic occultation light curve.

    Returns
    -------
    x_m : array
    intensity : array
    """

    # ─── Derived quantities ─────────────────────────────

    D_m = kbo.distance_au * AU_m
    R_m = kbo.radius_m
    b_m = kbo.impact_parameter_m

    # Star angular radius projected in the KBO, in meters
    star_radius_m = star.angular_radius_mas * mas_to_rad * D_m
    
    # Spatial grid
    x_m = np.linspace(-grid.x_max_m, grid.x_max_m, grid.n_x)

    # Define radial grid (must cover all possible r)
    r_max = np.sqrt((x_m.max() + star_radius_m)**2 + star_radius_m**2)
    r_grid_m = np.linspace(0, r_max, numerics.n_r_grid)

    # Wavelength grid
    lambdas_nm = np.linspace(
        bandpass.lam_min_nm,
        bandpass.lam_max_nm,
        bandpass.n_lambda
    )
    lambdas_m = lambdas_nm * nm_m

    # Spectral weights
    spec_w = planck_photon(lambdas_m, star.temperature_K)
    #filt_w = filter_transmission(
    #    lambdas_nm,
    #    bandpass.lam_min_nm,
    #    bandpass.lam_max_nm
    #)
    #weights = spec_w * filt_w
    
    # Fast loading of MAGIC transmission
    lam_qe, qe = load_response_file("optical_filter_MAGIC_QE.txt")
    qe_func = build_response_function(lam_qe, qe)
    response = qe_func
    if SII == True:
        lam_filt, filt = load_response_file("optical_filter_MAGIC_SII.txt")
        filt_func = build_response_function(lam_filt, filt)
        response = combine_responses(qe_func, filt_func)

    response_vals = response(lambdas_nm)
    weights = spec_w * response_vals

    #if bandpass.lam_min_nm != bandpass.lam_max_nm:
    #    weights /= weights.sum()
    weights /= weights.sum()

    # ─── Compute intensity ─────────────────────────────
    intensity_radial_total = np.zeros_like(r_grid_m)

    for lam, w in zip(lambdas_m, weights):
        I_r = fresnel_intensity_radial(
            r_grid_m,
            R_m,
            D_m,
            lam,
            n_int=numerics.n_int
        )
        intensity_radial_total += w * I_r

    if numerics.n_star_side > 1:
        intensity_total = apply_stellar_disk_2d(
            x_m,
            intensity_radial_total,
            r_grid_m,
            star_radius_m,
            b_m,
            numerics.n_star_side
        )
    return x_m, intensity_total

def simulate_poly_point(kbo, star, bandpass, grid, numerics):
    
    # KBO parameters
    D_m = kbo.distance_au * AU_m
    R_m = kbo.radius_m
    b_m = kbo.impact_parameter_m
    
    # Spatial grid
    x_m = np.linspace(-grid.x_max_m, grid.x_max_m, grid.n_x)

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
    return x_m, total

class OccultationEngine:
    """
    Precomputes static components and evaluates occultation light curves.
    """

    def __init__(self, star, bandpass, grid, numerics, response=None):

        self.star = star
        self.bandpass = bandpass
        self.numerics = numerics
        self.temperature_K = star.temperature_K
        self.response = response

        # --- wavelength grid ---
        self.lambdas_nm = np.linspace(
            bandpass.lam_min_nm,
            bandpass.lam_max_nm,
            bandpass.n_lambda
        )
        self.weights = self._compute_weights(response)

        # --- spatial grid ---
        self.x_m = np.linspace(-grid.x_max_m, grid.x_max_m, grid.n_x)

    def _compute_weights(self, response):
        spec_w = planck_photon(self.lambdas_nm * nm_m, self.temperature_K)

        if response is None:
            response_vals = filter_transmission(
                self.lambdas_nm,
                self.bandpass.lam_min_nm,
                self.bandpass.lam_max_nm
        )
        else:
            response_vals = response(self.lambdas_nm)

        weights = spec_w * response_vals
        return weights / weights.sum()

    def set_response(self, response=None):
        """
        Update instrument response and recompute spectral weights.

        Parameters
        ----------
        response : callable or None
            Function R(lambda_nm). If None, uses default bandpass.
        """
        self.response = response
        self.weights = self._compute_weights(response)

    def compute(self, kbo):
        """
        Compute occultation light curve for a given KBO.
        """
        
        # --- KBO parameters ---
        D_m = kbo.distance_au * AU_m
        R_m = kbo.radius_m
        b_m = kbo.impact_parameter_m

        # --- star projection ---
        r_star_m = self.star.angular_radius_mas * mas_to_rad * D_m

        if (self.star.angular_radius_mas < 0.0001):
            # For point source testing
            r_obs = np.sqrt(self.x_m**2 + b_m**2)
            intensity = np.zeros(len(self.x_m))
            
            for lam, w in zip(self.lambdas_nm * nm_m, self.weights):
                if w < 1e-12:
                    continue
                intensity += w * fresnel_intensity_radial(r_obs, R_m, D_m, lam, self.numerics.n_int)
            return self.x_m, intensity
        
        # --- radial grid ---
        r_max = np.sqrt(
            (self.x_m.max() + r_star_m)**2 +
            (r_star_m + b_m)**2
        )
        r_grid_m = np.linspace(0, r_max, self.numerics.n_r_grid)

        # --- polychromatic radial intensity ---
        intensity_radial = np.zeros_like(r_grid_m)

        for lam, w in zip(self.lambdas_nm * nm_m, self.weights):
            if w < 1e-12:
                # Skip this wavelength
                continue
            I_r = fresnel_intensity_radial(
                r_grid_m,
                R_m,
                D_m,
                lam,
                n_int=self.numerics.n_int
            )
            intensity_radial += w * I_r

        # --- 2D stellar convolution ---
        intensity = apply_stellar_disk_2d(
            self.x_m,
            intensity_radial,
            r_grid_m,
            r_star_m,
            b_m,
            self.numerics.n_star_side
        )

        return self.x_m, intensity

# ───────────────────────────────────────────────────────────
# Backward-compatible wrapper
# ───────────────────────────────────────────────────────────

def compute_lightcurve_test(kbo, star, bandpass, grid, numerics, response=None):
    """
    Legacy interface (kept for convenience).
    """

    engine = OccultationEngine(
        star,
        bandpass,
        grid,
        numerics,
        response=response,
    )

    return engine.compute(kbo)