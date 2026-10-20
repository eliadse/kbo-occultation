import numpy as np
import pytest

from kbo_occultation.simulation.kinematics import (ShadowKinematics, shadow_velocity_from_target,
                                        shadow_velocity_opposition)

astropy = pytest.importorskip("astropy")


def test_from_target_matches_analytic_at_opposition():
    from astropy.coordinates import SkyCoord
    from astropy.time import Time

    # A star close to the anti-solar point (opposition) on the ecliptic
    # in mid-January 2026.
    star = SkyCoord("07h40m00s +21d00m00s")
    t = Time("2026-01-15T22:00:00")
    kin = shadow_velocity_from_target(star, t, 43.0, include_kbo_orbit=True)
    analytic = shadow_velocity_opposition(43.0)
    assert abs(kin.v_rel_mps - analytic) / analytic < 0.05
    assert kin.direction_enu is None  # no location given


def test_earth_reflex_bounded_by_orbital_speed():
    from astropy.coordinates import SkyCoord
    from astropy.time import Time

    star = SkyCoord("02h47m00s +38d30m00s")  # off-ecliptic benchmark field
    t = Time("2026-01-15T22:00:00")
    kin = shadow_velocity_from_target(star, t, 43.0, include_kbo_orbit=False)
    assert 0.0 < kin.v_rel_mps < 30_500.0


def test_direction_enu_is_unit_vector():
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time

    loc = EarthLocation(lat=28.7619, lon=-17.89, height=2200)
    star = SkyCoord("07h40m00s +21d00m00s")
    t = Time("2026-01-15T22:00:00")
    kin = shadow_velocity_from_target(star, t, 43.0, location=loc)
    assert kin.direction_enu is not None
    assert np.isclose(np.linalg.norm(kin.direction_enu), 1.0)
