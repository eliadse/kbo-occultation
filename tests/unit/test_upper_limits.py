import numpy as np

from kbo_occultation.analysis.upper_limits import (AU_M, EfficiencyCurve, SurveyExposure,
                                          cumulative_density_upper_limit,
                                          effective_shadow_width_m, exposure_omega_deg2,
                                          fresnel_scale_m, poisson_upper_limit)
from kbo_occultation.simulation.kinematics import shadow_velocity_opposition


def test_poisson_upper_limit_values():
    assert np.isclose(poisson_upper_limit(0, 0.95), -np.log(0.05), rtol=1e-6)  # 2.9957
    assert np.isclose(poisson_upper_limit(0, 0.90), -np.log(0.10), rtol=1e-6)
    assert poisson_upper_limit(1, 0.95) > poisson_upper_limit(0, 0.95)


def test_fresnel_scale():
    # F = sqrt(lambda D / 2) ~ 1.27 km at 500 nm, 43 AU
    F = fresnel_scale_m(500e-9, 43.0)
    assert np.isclose(F, np.sqrt(500e-9 * 43.0 * AU_M / 2.0))
    assert 1200.0 < F < 1350.0


def test_effective_shadow_width_limits():
    F = fresnel_scale_m(500e-9, 43.0)
    # r -> 0: diffraction-dominated
    assert np.isclose(effective_shadow_width_m(0.0, 43.0, 500e-9), 2 * np.sqrt(3) * F)
    # star size adds 2 r_star
    assert np.isclose(
        effective_shadow_width_m(0.0, 43.0, 500e-9, star_radius_projected_m=100.0),
        2 * np.sqrt(3) * F + 200.0,
    )
    # r >> F: geometric limit 2r
    assert np.isclose(effective_shadow_width_m(1e6, 43.0, 500e-9), 2e6, rtol=1e-3)
    # monotonic in r
    r = np.linspace(0, 5000, 50)
    H = effective_shadow_width_m(r, 43.0, 500e-9)
    assert np.all(np.diff(H) >= 0)


def test_single_star_upper_limit_hand_computed():
    # One star, one radius, efficiency 1: N_UL = mu_UL / Omega.
    r0 = 500.0
    T = 3600.0
    v = 25_000.0
    D_au = 43.0
    lam = 500e-9

    eff = EfficiencyCurve(radii_m=np.array([100.0, 5000.0]), efficiency=np.array([1.0, 1.0]))
    exposure = SurveyExposure("HD-TEST", live_time_s=T, shadow_velocity_mps=v,
                              distance_au=D_au, efficiency=eff, wavelength_m=lam)

    D_m = D_au * AU_M
    H = effective_shadow_width_m(r0, D_au, lam)
    omega_deg2 = (H / D_m) * (v * T / D_m) * (180.0 / np.pi) ** 2
    expected = -np.log(0.05) / omega_deg2

    result = cumulative_density_upper_limit([exposure], np.array([r0]))
    assert np.isclose(result.n_upper_deg2[0], expected, rtol=1e-9)
    assert np.isclose(result.per_star_omega_deg2["HD-TEST"][0], omega_deg2, rtol=1e-12)
    assert np.isclose(exposure_omega_deg2(exposure, np.array([r0]))[0], omega_deg2)


def test_efficiency_zero_outside_measured_range():
    eff = EfficiencyCurve(radii_m=np.array([300.0, 1000.0]), efficiency=np.array([0.5, 1.0]))
    assert eff(100.0) == 0.0
    assert eff(2000.0) == 0.0
    assert np.isclose(eff(650.0), 0.75)

    exposure = SurveyExposure("S", 3600.0, 25_000.0, 43.0, eff)
    result = cumulative_density_upper_limit([exposure], np.array([100.0, 650.0]))
    assert np.isinf(result.n_upper_deg2[0])   # no sensitivity claimed there
    assert np.isfinite(result.n_upper_deg2[1])


def test_two_stars_combine():
    eff = EfficiencyCurve(radii_m=np.array([100.0, 5000.0]), efficiency=np.array([1.0, 1.0]))
    exp1 = SurveyExposure("A", 3600.0, 25_000.0, 43.0, eff)
    exp2 = SurveyExposure("B", 3600.0, 25_000.0, 43.0, eff)
    r = np.array([500.0])
    one = cumulative_density_upper_limit([exp1], r).n_upper_deg2[0]
    two = cumulative_density_upper_limit([exp1, exp2], r).n_upper_deg2[0]
    assert np.isclose(two, one / 2.0)


def test_opposition_velocity():
    v = shadow_velocity_opposition(43.0)
    assert 24_000.0 < v < 26_000.0
    # slower off opposition
    assert shadow_velocity_opposition(43.0, 40.0) < v
