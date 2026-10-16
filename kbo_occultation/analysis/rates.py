# kbo_occultation/analysis/rates.py
"""
KBO occultation *event-rate* equations (Nihei et al. 2007).

These are ported verbatim (physics-wise) from Víctor Pascual's exploratory
notebook ``tests/diffraction_pattern_expected_reduced.ipynb`` (where they were
named ``mu``/``velocity``/``surf_density``/``two_b``/``D_average``), cleaned up
and given consistent astropy-units handling.

They answer a different question from ``upper_limits.py``: given an assumed
sky surface density of KBOs, how *often* does a star get occulted, and hence
how long must you observe to expect one event. ``upper_limits.py`` is the
inverse (what density upper limit a null result implies for a given exposure).

The diffraction-aware shadow width H(r) (Nihei eq. 13) is NOT reimplemented
here; it is delegated to :func:`kbo_occultation.upper_limits.effective_shadow_width_m`
so there is a single source of truth. ``effective_shadow_width`` below is a
thin astropy-units wrapper around it.

References
----------
Nihei, T. C., et al. (2007), AJ 134, 1596 -- "Detectability of Occultations
of Stars by Objects in the Kuiper Belt and Oort Cloud". Equation numbers in
the docstrings refer to that paper.
"""

import numpy as np
from astropy import units as u
from astropy import constants as const

from .upper_limits import effective_shadow_width_m

# Size-distribution constants (Nihei et al. 2007, eq. 11 / 14).
Q_1 = 3.2e7        # normalisation of the cumulative surface density [deg^-2]
Q_L = 4.8          # large-object power-law slope


def _as_quantity(value, unit, name):
    """Accept a plain number (assumed already in `unit`) or a Quantity."""
    if isinstance(value, u.Quantity):
        if not value.unit.is_equivalent(unit):
            raise ValueError(f"Unexpected unit for {name}: {value.unit} (want {unit}).")
        return value
    return value * unit


def event_rate(shadow_width, velocity, surface_density, distance):
    """
    Occultation event rate for a single monitored star (Nihei eq. 8):

        mu = 2b * v * Sigma * (180 deg / pi / d)^2

    Parameters
    ----------
    shadow_width : Quantity [length]
        Effective shadow width 2b = H(r) (see :func:`effective_shadow_width`).
    velocity : Quantity [length/time]
        KBO shadow velocity across the observer (see :func:`kbo_velocity`, or
        ``kinematics.shadow_velocity_opposition``).
    surface_density : Quantity [1/deg^2]
        Sky surface density of KBOs of the relevant size (see
        :func:`surface_density`).
    distance : Quantity [length]
        Distance to the KBO.

    Returns
    -------
    Quantity [1/s]
        Expected occultation rate per star. The inverse, ``1 / mu``, is the
        mean observing time for one expected event (the y-axis of the
        "observing time vs KBO size" plot).
    """
    shadow_width = _as_quantity(shadow_width, u.m, "shadow width")
    velocity = _as_quantity(velocity, u.m / u.s, "velocity")
    surface_density = _as_quantity(surface_density, 1 / u.deg**2, "surface density")
    distance = _as_quantity(distance, u.m, "KBO distance")

    rate = shadow_width * velocity * surface_density * (180 * u.deg / np.pi / distance) ** 2
    return rate.to(1 / u.s)


def kbo_velocity(r_earth, r_kbo, elongation):
    """
    Transverse (shadow) velocity of a KBO on a circular orbit, as projected
    onto the observer, as a function of solar elongation (Nihei eq. 10):

        v = sqrt(G M_sun / r_earth) *
            [ sqrt( (r_e/r_k) * (1 - (r_e/r_k)^2 sin^2 e) ) + cos(e) ]

    At opposition (elongation = 180 deg) this is the maximum retrograde
    (Earth-reflex-dominated) shadow speed, ~25 km/s near 40 AU.

    Parameters
    ----------
    r_earth, r_kbo : Quantity [length]
        Orbital radii of Earth and the KBO.
    elongation : Quantity [angle]
        Solar elongation of the field.

    Returns
    -------
    Quantity [m/s]
    """
    r_earth = _as_quantity(r_earth, u.m, "Earth orbital radius")
    r_kbo = _as_quantity(r_kbo, u.m, "KBO orbital radius")
    elongation = _as_quantity(elongation, u.deg, "elongation")

    prefactor = (const.G * const.M_sun / r_earth) ** 0.5
    ratio = r_earth / r_kbo
    factor = (ratio * (1 - ratio**2 * np.sin(elongation) ** 2)) ** 0.5 + np.cos(elongation)
    return (prefactor * factor).to(u.m / u.s)


def surface_density(diameter_break, diameter, q_s):
    """
    Cumulative sky surface density of KBOs at diameter `diameter`
    (Nihei et al. 2007, eq. 11):

        Sigma = Q_1 * (D_break[km])^(q_s - q_L) * (D[km])^(1 - q_s)   [deg^-2]

    with Q_1 = 3.2e7, q_L = 4.8 (module constants).

    Parameters
    ----------
    diameter_break : Quantity [length]
        Break diameter D_k (Nihei uses 80 km).
    diameter : Quantity [length]
        Diameter(s) at which to evaluate the density (D_0 in the notebook).
    q_s : float
        Small-object power-law slope of the size distribution.

    Returns
    -------
    Quantity [1/deg^2]
    """
    q_s = float(q_s.value if isinstance(q_s, u.Quantity) else q_s)
    diameter_break = _as_quantity(diameter_break, u.m, "break diameter")
    diameter = _as_quantity(diameter, u.m, "diameter")

    sigma = (Q_1
             * diameter_break.to(u.km).value ** (q_s - Q_L)
             * diameter.to(u.km).value ** (1 - q_s))
    return sigma / u.deg**2


def mean_diameter(diameter_break, diameter, q_s):
    """
    Effective (flux-weighted mean) KBO diameter for a broken power-law size
    distribution (Nihei et al. 2007, eq. 14). This is the single "effective
    size" fed to :func:`effective_shadow_width` for the smooth analytic
    ("parameterization") rate curve.

    Parameters
    ----------
    diameter_break : Quantity [length]
        Break diameter D_k.
    diameter : Quantity [length]
        Reference diameter D_0.
    q_s : float
        Small-object slope.

    Returns
    -------
    Quantity [length]
    """
    q_s = float(q_s.value if isinstance(q_s, u.Quantity) else q_s)
    diameter_break = _as_quantity(diameter_break, u.m, "break diameter")
    diameter = _as_quantity(diameter, u.m, "diameter")

    if q_s == 2:
        factor = (1 - Q_L) / (2 - Q_L) + np.log(diameter_break / diameter)
    else:
        factor = (((1 - Q_L) / (2 - Q_L) - (1 - q_s) / (2 - q_s))
                  * (diameter / diameter_break) ** (q_s - 2)
                  + (1 - q_s) / (2 - q_s))
    return (factor * diameter).to(u.m)


def effective_shadow_width(wavelength, distance, kbo_radius, star_radius):
    """
    Diffraction-aware effective shadow width 2b = H(r), Nihei eq. 13 --
    astropy-units wrapper around
    :func:`kbo_occultation.upper_limits.effective_shadow_width_m`:

        H(r) = 2 * ((sqrt(3) F)^1.5 + r^1.5)^(2/3) + 2 r_star,
        F = sqrt(lambda d / 2)

    (The notebook wrote the equivalent ``((sqrt(6 lambda d))^1.5 + r^1.5)^(2/3)
    + 2 r_star``; sqrt(6 lambda d) == 2 sqrt(3) F, so the two forms are
    algebraically identical.)

    Parameters
    ----------
    wavelength : Quantity [length]
        Effective observing wavelength.
    distance : Quantity [length]
        Distance to the KBO.
    kbo_radius : Quantity [length]
        KBO radius r (scalar or array).
    star_radius : Quantity [length]
        Projected stellar radius at the KBO distance
        (= tan(theta_star) * distance).

    Returns
    -------
    Quantity [m]
    """
    wavelength = _as_quantity(wavelength, u.m, "wavelength")
    distance = _as_quantity(distance, u.m, "KBO distance")
    kbo_radius = _as_quantity(kbo_radius, u.m, "KBO radius")
    star_radius = _as_quantity(star_radius, u.m, "stellar radius")

    distance_au = distance.to(u.au).value
    width_m = effective_shadow_width_m(
        kbo_radius.to(u.m).value,
        distance_au,
        wavelength.to(u.m).value,
        star_radius.to(u.m).value,
    )
    return width_m * u.m
