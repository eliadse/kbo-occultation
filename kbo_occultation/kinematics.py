# kbo_occultation/kinematics.py
"""
Shadow-velocity helpers.

The KBO shadow sweeps across the observer at the relative transverse
velocity between the (hypothetical) occulter and the Earth. Two
contexts in this package need it, with different conventions
(deliberately -- see the campaign design):

- Coincidence timing between telescopes: use the *actual* kinematics of
  the observation. For off-ecliptic benchmark stars there is no real
  KBO, so the well-defined choice is the Earth's reflex motion alone
  (`shadow_velocity_from_target(..., include_kbo_orbit=False)`).
- Injection stretching and exposure bookkeeping: use the velocity the
  benchmark run is meant to *emulate* -- typically the opposition value
  (`shadow_velocity_opposition`), so the measured efficiency curve
  carries over directly to future on-ecliptic observations.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np

EARTH_ORBITAL_SPEED_MPS = 29_784.0  # mean heliocentric speed of the Earth


@dataclass
class ShadowKinematics:
    """
    v_rel_mps : relative shadow speed across the observer (m/s).
    direction_enu : optional unit 3-vector (east, north, up) of the
        shadow's drift direction in the local frame. None means unknown,
        which treats the telescope baselines isotropically (conservative).
    """
    v_rel_mps: float
    direction_enu: Optional[np.ndarray] = None


def shadow_velocity_opposition(distance_au: float, opposition_angle_deg: float = 0.0) -> float:
    """
    Retrograde shadow speed for a KBO on a circular ecliptic orbit at
    ``distance_au``, observed at ``opposition_angle_deg`` from
    opposition (Nihei et al. 2007, the convention used by TAOS):

        v_rel = v_E * (cos(eps) - a**-0.5 * sqrt(1 - a**-2 * sin(eps)**2))

    At opposition (eps = 0) and a = 43 AU this gives ~25.2 km/s.
    """
    a = float(distance_au)
    eps = np.deg2rad(opposition_angle_deg)
    return EARTH_ORBITAL_SPEED_MPS * (
        np.cos(eps) - a**-0.5 * np.sqrt(1.0 - a**-2 * np.sin(eps)**2)
    )


def shadow_velocity_from_target(star_coord, obs_time, distance_au: float,
                                include_kbo_orbit: bool = True,
                                location=None) -> ShadowKinematics:
    """
    Shadow kinematics computed from the actual observing geometry, valid
    for any star (including off-ecliptic benchmark fields).

    The shadow of a star (at infinity) cast by an occulter moves across
    the Earth with the transverse (perpendicular-to-line-of-sight)
    relative velocity v_perp(occulter) - v_perp(Earth).

    Parameters
    ----------
    star_coord : astropy.coordinates.SkyCoord
        Target star (ICRS).
    obs_time : astropy.time.Time
        Time of observation.
    distance_au : float
        Distance at which an occulter is placed on the line of sight.
    include_kbo_orbit : bool
        If True, subtract the transverse velocity of a hypothetical occulter 
        on a prograde circular heliocentric orbit through that point (orbit 
        plane chosen so its velocity is parallel to the ecliptic plane -- the 
        natural extension of the on-ecliptic convention). If False, use the 
        Earth's reflex motion alone -- the well-defined choice for benchmark
        stars where no occulter exists (used for coincidence timing).
    location : astropy.coordinates.EarthLocation, optional
        If given, the shadow drift direction is also expressed in the local ENU
        frame (for projecting telescope baselines); otherwise ``direction_enu``
        is None and baselines must be treated isotropically.

    Returns
    -------
    ShadowKinematics
    """
    from astropy import units as u
    from astropy.coordinates import (AltAz, SkyCoord, get_body_barycentric_posvel,
                                     CartesianRepresentation, BarycentricMeanEcliptic)

    # Line of sight (unit vector, ICRS cartesian)
    los = star_coord.icrs.cartesian
    n_hat = (los / los.norm()).get_xyz().value

    earth_pos, earth_vel = get_body_barycentric_posvel("earth", obs_time)
    v_earth = earth_vel.get_xyz().to_value(u.m / u.s)

    def perp(v):
        return v - np.dot(v, n_hat) * n_hat

    v_rel_vec = -perp(v_earth)

    if include_kbo_orbit:
        # Hypothetical occulter on the LOS at distance_au from Earth.
        AU_m = u.au.to(u.m)
        r_obj = earth_pos.get_xyz().to_value(u.m) + distance_au * AU_m * n_hat
        r_obj_norm = np.linalg.norm(r_obj)
        # Circular-orbit speed at that heliocentric radius.
        GM_sun = 1.32712440018e20  # m^3 / s^2
        v_circ = np.sqrt(GM_sun / r_obj_norm)
        # Prograde direction parallel to the ecliptic plane: z_ecl x r_hat.
        # (Barycentric ecliptic frame: pure rotation w.r.t. ICRS, so a unit
        # direction transforms without needing a distance.)
        z_ecl_icrs = SkyCoord(
            CartesianRepresentation(0.0, 0.0, 1.0),
            frame=BarycentricMeanEcliptic(equinox="J2000"),
        ).icrs.cartesian.get_xyz().value
        z_ecl_icrs /= np.linalg.norm(z_ecl_icrs)
        v_dir = np.cross(z_ecl_icrs, r_obj / r_obj_norm)
        v_dir /= np.linalg.norm(v_dir)
        v_rel_vec = perp(v_circ * v_dir) - perp(v_earth)

    v_rel = float(np.linalg.norm(v_rel_vec))

    direction_enu = None
    if location is not None and v_rel > 0:
        # Express the drift direction in the local ENU frame by treating it as a 
        # direction on the sky (directions transform like positions at infinity).
        d_icrs = SkyCoord(CartesianRepresentation(*(v_rel_vec / v_rel)), frame="icrs")
        altaz = d_icrs.transform_to(AltAz(obstime=obs_time, location=location))
        alt = altaz.alt.rad
        az = altaz.az.rad
        direction_enu = np.array([
            np.cos(alt) * np.sin(az),   # east
            np.cos(alt) * np.cos(az),   # north
            np.sin(alt),                # up
        ])

    return ShadowKinematics(v_rel_mps=v_rel, direction_enu=direction_enu)
