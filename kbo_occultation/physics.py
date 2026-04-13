# kbo_occultation/physics.py

import numpy as np
from scipy.special import j0


# ─── Physical constants ─────────────────────────────────────

AU_m = 1.49598e11
km_m = 1.0e3
nm_m = 1.0e-9
mas_to_rad = np.pi / (180.0 * 3600.0 * 1.0e3)


# ─── Spectral functions ─────────────────────────────────────

def planck_photon(lam_m, T_K):
    """Planck photon-count spectrum."""
    h, c, kB = 6.626e-34, 2.998e8, 1.381e-23
    x = h * c / (lam_m * kB * T_K)
    x = np.clip(x, 1e-30, 700)
    return lam_m**(-4) / (np.exp(x) - 1.0)


def filter_transmission(lam_nm, lam_min_nm, lam_max_nm):
    """Dummy filter, a simple trapezoidal bandpass."""
    lam = np.asarray(lam_nm, dtype=float)
    edge = 5.0

    T = np.clip((lam - lam_min_nm) / edge, 0.0, 1.0)
    T *= np.clip((lam_max_nm - lam) / edge, 0.0, 1.0)
    return T


# ─── Fresnel diffraction ────────────────────────────────────

def fresnel_intensity_radial(r_arr_m, R_m, D_m, lam_m, n_int=800):
    """
    Fresnel diffraction pattern of a circular opaque disk.
    Uses Lommel–Hankel integral form.
    """

    # Fresnel scale
    F = np.sqrt(lam_m * D_m / 2.0)

    r_F = r_arr_m / F
    R_F = R_m / F

    # Integration grid
    t = np.linspace(0, R_F, n_int)
    dt = t[1] - t[0]

    intensity = np.zeros_like(r_arr_m)

    for i, rF in enumerate(r_F):
        integrand = t * j0(np.pi * t * rF) * np.exp(1j * np.pi * t**2 / 2)
        u_disk = -1j * np.exp(1j * np.pi * rF**2 / 2) * 2 * np.pi * np.sum(integrand) * dt
        U = 1 - u_disk
        intensity[i] = np.abs(U)**2

    return intensity
