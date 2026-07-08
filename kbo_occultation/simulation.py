# kbo_occultation/simulation.py

import numpy as np

from .physics import (planck_photon, filter_transmission, fresnel_intensity_radial, AU_m, km_m, nm_m, mas_to_rad, fresnel_point_intensity)
from .instruments import *

def simulate_poly_point(kbo, star, bandpass, grid, numerics, weight_rel_tol=1e-6):
    
    # KBO parameters
    D_m = kbo.distance_au * AU_m
    R_m = kbo.radius_m
    b_m = kbo.impact_parameter_m
    
    # Spatial grid
    x_m = np.linspace(-grid.x_max_m, grid.x_max_m, grid.n_x)

    # Wavelength grid
    lambdas_nm = np.linspace(bandpass.lam_min_nm, bandpass.lam_max_nm, bandpass.n_lambda)
    lambdas_m = lambdas_nm * nm_m

    # Spectral weights
    spec_w = planck_photon(lambdas_m, star.temperature_K)
    filt_w = filter_transmission(lambdas_nm, bandpass.lam_min_nm, bandpass.lam_max_nm)
    weights = spec_w * filt_w
    weights /= weights.sum()

    N_int = numerics.n_int

    """ Monochromatic and polychromatic point source"""
    r_obs = np.sqrt(x_m**2 + b_m**2)
    total = np.zeros(len(x_m))
    weight_floor = weight_rel_tol * weights.max()
    for lam_m, w in zip(lambdas_m, weights):
        if w < weight_floor:
            continue
        #total += w * fresnel_intensity_radial(r_obs, R_m, D_m, lam_m, N_int)
        F = np.sqrt(lam_m * D_m / 2.0)

        r = r_obs / F
        rho = R_m / F
        total += w * fresnel_point_intensity(r, rho, N_int)
    return x_m, total

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
        shifted = np.interp(x_m + shift, x_m, intensity, left=1.0, right=1.0)
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

def compute_lightcurve_old(kbo, star, bandpass, grid, numerics, SII=False):
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
    lambdas_nm = np.linspace(bandpass.lam_min_nm, bandpass.lam_max_nm, bandpass.n_lambda)
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
        F = np.sqrt(lam * D_m / 2.0)

        r = r_grid_m / F
        rho = R_m / F
        I_r = fresnel_point_intensity(r, rho, n_int=numerics.n_int)
        #I_r = fresnel_intensity_radial(r_grid_m, R_m, D_m, lam, n_int=numerics.n_int)
        intensity_radial_total += w * I_r

    if numerics.n_star_side > 1:
        intensity_total = apply_stellar_disk_2d(x_m, intensity_radial_total, r_grid_m, star_radius_m, b_m, numerics.n_star_side)
    return x_m, intensity_total

class OccultationEngine:
    """
    Precomputes static components and evaluates occultation light curves.

    Design note on caching
    -----------------------
    The per-wavelength Fresnel radial intensity profile I(r) depends only
    on the KBO radius and distance (R_m, D_m) -- NOT on impact parameter
    or stellar angular size, which only enter later, in the stellar-disk
    convolution step. To exploit that, this engine evaluates the profile
    on a *fixed* radial grid (``self.r_grid_m``), sized once at
    construction time to safely cover every impact parameter / stellar
    size you intend to pass to ``compute()``. That fixed grid is what
    makes it safe to cache profiles by ``(R_m, D_m)`` alone and reuse
    them across a whole parameter sweep (see ``sweep.py``).

    If you call ``compute()`` with a star size or impact parameter that
    needs a larger radial grid than the engine was built for, the grid is
    transparently extended (and the profile cache is invalidated, since
    it was keyed for the old grid) -- with a warning, since this silently
    costs you the reuse benefit for the rest of the sweep. Pass a generous
    ``r_max_m`` up front (or use ``sweep.run_parameter_sweep``, which
    sizes it automatically) to avoid this.
    """

    def __init__(self, star, bandpass, grid, numerics, response=None, r_max_m=None, r_max_margin=1.2):

        self.star = star
        self.bandpass = bandpass
        self.numerics = numerics
        self.temperature_K = star.temperature_K
        self.response = response

        # --- wavelength grid ---
        self.lambdas_nm = np.linspace(bandpass.lam_min_nm, bandpass.lam_max_nm, bandpass.n_lambda)
        self.weights = self._compute_weights(response)

        # --- spatial grid ---
        self.x_m = np.linspace(-grid.x_max_m, grid.x_max_m, grid.n_x)

        # --- fixed radial grid (see class docstring) ---
        if r_max_m is None:
            # Covers a point source at zero impact parameter only. Widen
            # this (via r_max_m) before sweeping over star size or impact
            # parameter, or let sweep.run_parameter_sweep size it for you.
            r_max_m = r_max_margin * self.x_m.max()
        self.r_max_m = r_max_m
        self.r_grid_m = np.linspace(0, self.r_max_m, self.numerics.n_r_grid)

        # --- cached fresnel profiles, keyed by (R_m, D_m) only, since
        #     r_grid_m is now fixed for the life of this engine ---
        self._fresnel_cache = {}

    def _compute_weights(self, response):
        spec_w = planck_photon(self.lambdas_nm * nm_m, self.temperature_K)

        if response is None:
            response_vals = filter_transmission(self.lambdas_nm, self.bandpass.lam_min_nm, self.bandpass.lam_max_nm)
        else:
            response_vals = response(self.lambdas_nm)

        weights = spec_w * response_vals
        return weights / weights.sum()

    def _ensure_radial_grid_covers(self, r_needed_m):
        """
        Grow the fixed radial grid (and invalidate the profile cache) if a
        request needs more range than the engine was built for.
        """
        if r_needed_m <= self.r_max_m:
            return
        import warnings
        warnings.warn(
            f"OccultationEngine radial grid (r_max_m={self.r_max_m:.3g}) was "
            f"too small for this call (needed {r_needed_m:.3g} m); extending "
            "it and clearing the Fresnel-profile cache. Pass a larger "
            "r_max_m at construction time to avoid recomputing profiles "
            "mid-sweep.",
            stacklevel=3,
        )
        self.r_max_m = r_needed_m * 1.05
        self.r_grid_m = np.linspace(0, self.r_max_m, self.numerics.n_r_grid)
        self._fresnel_cache = {}

    def _get_fresnel_profiles(self, R_m, D_m):
        """
        Returns cached per-wavelength Fresnel profiles I(r), evaluated on
        self.r_grid_m, for all wavelengths in this engine's bandpass.
        """

        key = (R_m, D_m)

        if key in self._fresnel_cache:
            return self._fresnel_cache[key]

        profiles = []

        for lam in self.lambdas_nm * nm_m:
            F = np.sqrt(lam * D_m / 2.0)
            r = self.r_grid_m / F
            rho = R_m / F
            I_r = fresnel_point_intensity(r, rho, self.numerics.n_int)

            profiles.append(I_r)

        profiles = np.array(profiles)  # shape (n_lambda, n_r)

        self._fresnel_cache[key] = profiles

        return profiles

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

    def get_response(self):
        """
        Get the weights that the engine is currently using.
        """
        return self.lambdas_nm, self.weights

    def compute(self, kbo, star_angular_radius_mas=None, weight_rel_tol=1e-6):
        """
        Compute occultation light curve for a given KBO.

        Parameters
        ----------
        kbo : KBOConfig
        star_angular_radius_mas : float, optional
            Override the engine's star angular size for this call only,
            without rebuilding the engine or losing the Fresnel-profile
            cache. This is what lets a parameter sweep vary stellar size
            (and impact parameter, via kbo.impact_parameter_m) cheaply
            across many calls to the same engine. Defaults to
            self.star.angular_radius_mas.
        weight_rel_tol : float, optional
            Wavelengths whose spectral weight is below this fraction of
            the bandpass's maximum weight are skipped in the point-source
            branch. Relative (not absolute) so this stays sensible for
            narrow filters, where every weight can be individually small.
            Default 1e-6.

        Returns
        -------
        x_m, intensity : ndarray, ndarray
        """

        # --- KBO parameters ---
        D_m = kbo.distance_au * AU_m
        R_m = kbo.radius_m
        b_m = kbo.impact_parameter_m

        # --- star projection ---
        star_mas = (self.star.angular_radius_mas if star_angular_radius_mas is None
                    else star_angular_radius_mas)
        r_star_m = star_mas * mas_to_rad * D_m

        if star_mas < 0.0001:
            # Point source: evaluate directly at the observer offsets, no
            # radial grid or stellar-disk convolution needed.
            r_obs = np.sqrt(self.x_m**2 + b_m**2)
            intensity = np.zeros(len(self.x_m))

            weight_floor = weight_rel_tol * self.weights.max()
            for lam, w in zip(self.lambdas_nm * nm_m, self.weights):
                if w < weight_floor:
                    continue

                F = np.sqrt(lam * D_m / 2.0)
                r = r_obs / F
                rho = R_m / F
                intensity += w * fresnel_point_intensity(r, rho, self.numerics.n_int)
            return self.x_m, intensity

        # --- radial grid: extend the fixed grid if this call needs more
        #     range than the engine was built for (see class docstring) ---
        r_needed = np.sqrt((self.x_m.max() + r_star_m)**2 + (r_star_m + b_m)**2)
        self._ensure_radial_grid_covers(r_needed)

        # --- polychromatic radial intensity, on the shared fixed grid ---
        profiles = self._get_fresnel_profiles(R_m, D_m)

        # Apply weights from transmission/filters etc. Weighted sum over wavelength axis
        intensity_radial = np.tensordot(self.weights, profiles, axes=(0, 0))

        # --- 2D stellar convolution ---
        intensity = apply_stellar_disk_2d(
            self.x_m, intensity_radial, self.r_grid_m, r_star_m, b_m, self.numerics.n_star_side
        )

        return self.x_m, intensity

# ───────────────────────────────────────────────────────────
# Backward-compatible wrapper
# ───────────────────────────────────────────────────────────

def compute_lightcurve(kbo, star, bandpass, grid, numerics, response=None, r_max_m=None):
    """
    Legacy interface (kept for convenience): builds a fresh, one-shot
    engine sized exactly for this call.

    For repeated calls -- e.g. a parameter sweep over impact parameter,
    stellar size, or KBO radius/distance -- build an OccultationEngine
    yourself and reuse it (or use sweep.run_parameter_sweep), so the
    Fresnel-profile cache actually pays off. See OccultationEngine's
    docstring for why a shared engine needs a shared radial grid.
    """
    if r_max_m is None:
        D_m = kbo.distance_au * AU_m
        r_star_m = star.angular_radius_mas * mas_to_rad * D_m
        b_m = kbo.impact_parameter_m
        r_max_m = np.sqrt((grid.x_max_m + r_star_m)**2 + (r_star_m + b_m)**2)

    engine = OccultationEngine(star, bandpass, grid, numerics, response=response, r_max_m=r_max_m)

    return engine.compute(kbo)