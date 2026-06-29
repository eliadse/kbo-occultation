# kbo_occultation/physics.py

import numpy as np
from scipy.special import j0
from scipy.special import jn


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


def lommel_U(n, mu, nu, iterations=40):
    """
    Evaluate the Lommel function U_n(mu, nu) used in Fresnel diffraction
    by a circular opaque disk.

    The Lommel functions are a classical series solution to the Fresnel
    diffraction integral for a circular aperture/obstacle. They appear
    naturally when the integral is expanded in a power series of Bessel
    functions.

    Definition:
        U_n(mu, nu) = sum_{k=0}^{inf} (-1)^k * (mu/nu)^(n+2k) * J_{n+2k}(pi*mu*nu)

    where J_m is the Bessel function of the first kind of order m.

    IMPORTANT: This series converges only when mu/nu < 1. The calling
    function (intensity) is responsible for passing arguments in the
    correct order depending on the geometric region (see intensity() below).

    Parameters
    ----------
    n : int
        Order of the Lommel function (0, 1, or 2 in practice).
    mu : float or array_like
        The smaller of the two radii (in Fresnel units), so that mu/nu < 1.
    nu : float or array_like
        The larger of the two radii (in Fresnel units).
    iterations : int, optional
        Number of terms to include in the series. Default is 40.
        More terms are needed when mu/nu is close to 1 (near the shadow edge).

    Returns
    -------
    float or ndarray
        Value of U_n(mu, nu).

    References
    ----------
    Lacki, B. C. (2014) - analytic Fresnel diffraction for a circular disk.
    Born & Wolf, "Principles of Optics", Chapter 8.
    """
    return np.sum(
        [
            np.power(-1, k)
            * np.power(mu / nu, n + (2 * k))
            * jn(n + (2 * k), np.pi * mu * nu)
            for k in np.arange(0, iterations, 1)
        ],
        axis=0,
    )


def _lommel_U0_edge(rho):
    """
    Evaluate U_0 at the geometric shadow edge (r ≈ rho) using the
    Fresnel-integral form.

    At the shadow edge the standard Lommel series (lommel_U) is poorly
    conditioned because mu/nu → 1, making the series converge very slowly.
    This closed-form expression avoids that problem.

    Definition:
        U_0(rho) = 0.5 * [cos(pi * rho^2) + J_0(pi * rho^2)]

    Parameters
    ----------
    rho : float or array_like
        Disk radius in Fresnel units (= r at the edge, since r ≈ rho).

    Returns
    -------
    float or ndarray
        Value of U_0 at the shadow edge.
    """
    return 0.5 * (np.cos(np.pi * rho**2) + jn(0, np.pi * rho**2))


def _lommel_U1_edge(rho):
    """
    Evaluate U_1 at the geometric shadow edge (r ≈ rho) using the
    Fresnel-integral form.

    Companion to _lommel_U0_edge; together they give the intensity at
    the edge as U_0^2 + U_1^2.

    Definition:
        U_1(rho) = 0.5 * sin(pi * rho^2)

    Parameters
    ----------
    rho : float or array_like
        Disk radius in Fresnel units.

    Returns
    -------
    float or ndarray
        Value of U_1 at the shadow edge.
    """
    return 0.5 * np.sin(np.pi * rho**2)


def fresnel_point_intensity(r, rho, iterations=40):
    """
    Compute the Fresnel diffraction intensity of a monochromatic point source
    behind a circular opaque disk, normalised to the unobscured flux.

    Uses the Lommel-function series solution (Lacki 2014), which is the
    analytic equivalent of the Huygens-Fresnel integral. All lengths are
    expressed in Fresnel units:

        r_F  =  r_physical / sqrt(lambda * D / 2)

    where D is the observer-to-screen distance and lambda the wavelength.

    The solution is split into three geometric regions because the Lommel
    series U_n(mu, nu) only converges when mu/nu < 1.  In each region the
    arguments are ordered so that the smaller quantity is always passed as
    mu and the larger as nu:

    Region 1 — inside the geometric shadow (r < rho):
        I = U_n(0, r, rho)^2 + U_n(1, r, rho)^2
        (mu = r, nu = rho, so mu/nu < 1 ✓)

    Region 2 — outside the geometric shadow (r > rho):
        I = 1 + U_n(1, rho, r)^2 + U_n(2, rho, r)^2
              - 2*U_n(1, rho, r)*sin(pi/2 * (rho^2 + r^2))
              + 2*U_n(2, rho, r)*cos(pi/2 * (rho^2 + r^2))
        (mu = rho, nu = r, so mu/nu < 1 ✓)
        The leading "1" is the unobstructed wave; the remaining terms
        are the diffraction correction.

    Region 3 — at the geometric shadow edge (r ≈ rho):
        Lommel series converges too slowly; use the closed-form edge
        expressions _lommel_U0_edge and _lommel_U1_edge instead.

    Parameters
    ----------
    r : array_like
        Radial offset of the observer from the line of sight, in Fresnel
        units. Negative values are treated as positive (radial symmetry).
    rho : float
        Radius of the opaque disk in Fresnel units.
    iterations : int, optional
        Number of terms in the Lommel series. Default 40 is sufficient for
        |r - rho| / rho > ~0.001. Increase if you need accuracy very close
        to the shadow edge.

    Returns
    -------
    ndarray
        Intensity normalised to the unobscured source flux (I/I_0).
        Values near r=0 approach 1.0 (the Poisson/Arago bright spot).

    Notes
    -----
    * The rtol=1e-3 threshold for the edge region is a practical choice;
      tighten it (smaller rtol) if you observe discontinuities near r=rho.
    * For large r or rho (many Fresnel units) you may need more iterations
      because higher-order Bessel functions decay more slowly.

    References
    ----------
    Lacki, B. C. (2014), "On the use of Cherenkov Telescopes for outer Solar 
    system body occultations".
    Born & Wolf, "Principles of Optics", 7th ed., §8.7.
    """
    r = np.abs(np.asarray(r, dtype=float))
    intensity = np.zeros_like(r)

    # Boolean masks for the three regions
    at_edge    = np.isclose(r, rho, rtol=1e-03)
    in_shadow  = np.logical_and(r < rho, ~at_edge)
    out_shadow = np.logical_and(r > rho, ~at_edge)

    # --- Region 3: shadow edge ---
    # Use the special closed-form U_0, U_1 to avoid slow series convergence.
    if np.any(at_edge):
        r_e = r[at_edge]
        intensity[at_edge] = (_lommel_U0_edge(r_e) ** 2 + _lommel_U1_edge(r_e) ** 2)

    # --- Region 1: inside the geometric shadow (r < rho) ---
    # Lommel series with mu=r, nu=rho ensures mu/nu < 1.
    if np.any(in_shadow):
        r_i = r[in_shadow]
        intensity[in_shadow] = (
            lommel_U(0, r_i, rho, iterations) ** 2
            + lommel_U(1, r_i, rho, iterations) ** 2
        )

    # --- Region 2: outside the geometric shadow (r > rho) ---
    # Lommel series with mu=rho, nu=r ensures mu/nu < 1.
    # The formula comes from decomposing the total field into the
    # unobstructed wave (=1) plus the diffracted contribution from
    # the disk edge, which interferes constructively/destructively
    # depending on the phase term (rho^2 + r^2).
    if np.any(out_shadow):
        r_o = r[out_shadow]
        U1 = lommel_U(1, rho, r_o, iterations)
        U2 = lommel_U(2, rho, r_o, iterations)
        phase = (np.pi / 2) * (rho**2 + r_o**2)

        intensity[out_shadow] = (
            1
            + U1**2
            + U2**2
            - 2 * U1 * np.sin(phase)
            + 2 * U2 * np.cos(phase)
        )

    return intensity
