import numpy as np
import pytest
from astropy import units as u

from kbo_occultation.analysis import rates
from kbo_occultation.analysis.upper_limits import EfficiencyCurve

V = 25e3 * u.m / u.s
D = (39 * u.au).to(u.m)


def test_observing_time_per_size_is_the_inverse_of_event_rate():
    r = np.array([100.0, 300.0])
    w = np.array([3000.0, 3200.0])
    t = rates.observing_time_per_size(r, w, V, D, 80 * u.km, 4.0)
    for ri, wi, ti in zip(r, w, t):
        mu = rates.event_rate(wi * u.m, V, rates.surface_density(80 * u.km, 2 * ri * u.m, 4.0), D)
        assert ti.to(u.s).value == pytest.approx((1 / mu).to(u.s).value, rel=1e-12)
    assert np.isinf(rates.observing_time_per_size([100.0], [0.0], V, D).value[0])


def test_detected_event_rate_of_a_constant_width_is_the_cumulative_density():
    r = np.geomspace(20.0, 1000.0, 30)
    w = np.full_like(r, 3000.0)
    mu = rates.detected_event_rate(r, w, V, D, 80 * u.km, 4.0, extrapolate=False)
    sigma = (rates.surface_density(80 * u.km, 2 * r[0] * u.m, 4.0)
             - rates.surface_density(80 * u.km, 2 * r[-1] * u.m, 4.0))
    assert mu.to(1 / u.s).value == pytest.approx(rates.event_rate(3000 * u.m, V, sigma, D).to(1 / u.s).value, rel=1e-3)


def test_detected_event_rate_ignores_undetectable_sizes_and_adds_the_large_ones():
    r = np.geomspace(20.0, 1000.0, 30)
    w = np.where(r >= 100.0, 3000.0, 0.0)
    all_sizes = rates.detected_event_rate(r, np.full_like(r, 3000.0), V, D, extrapolate=False)
    big_only = rates.detected_event_rate(r, w, V, D, extrapolate=False)
    assert 0.0 < big_only.value < all_sizes.value
    cut = rates.detected_event_rate(r[r >= 100.0], w[r >= 100.0], V, D, extrapolate=False)
    assert big_only.value == pytest.approx(cut.value, rel=0.02)
    with_tail = rates.detected_event_rate(r, 2.0 * r + 2000.0, V, D, extrapolate=True)
    without = rates.detected_event_rate(r, 2.0 * r + 2000.0, V, D, extrapolate=False)
    assert without.value < with_tail.value < 1.05 * without.value          # the large objects are a few per cent


def test_efficiency_curve_clip_flag():
    r, e = np.array([10.0, 100.0]), np.array([0.5, 1.5])
    assert EfficiencyCurve(r, e)(100.0) == 1.0
    assert EfficiencyCurve(r, e, clip=False)(100.0) == pytest.approx(1.5)
    assert EfficiencyCurve(r, e, clip=False)(1000.0) == 0.0
