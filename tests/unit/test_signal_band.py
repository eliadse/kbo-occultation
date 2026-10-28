import numpy as np
import pytest

from kbo_occultation import BandpassConfig, GridConfig, NumericalConfig
from kbo_occultation.analysis.bands import signal_band
from kbo_occultation.simulation.detectability import (
    cumulative_quantile_frequency, spatial_power_spectrum)

BAND = BandpassConfig(400, 500, 11)
NUM = NumericalConfig(n_r_grid=3000, n_star_side=48)


def _band(star=0.02, v=25e3, grid=GridConfig(20000, 4000), **kw):
    return signal_band(star, [100], BAND, 40.0, v, grid=grid, numerics=NUM,
                       temperature_K=10000.0, **kw)[0]


def test_frequency_scales_linearly_with_velocity():
    a, b = _band(v=10e3), _band(v=30e3)
    for key in ("f05_hz", "f50_hz", "f90_hz", "f99_hz"):
        assert b[key] / a[key] == pytest.approx(3.0, rel=1e-6)


def test_stable_under_grid_refinement():
    coarse = _band(grid=GridConfig(20000, 4000))
    fine = _band(grid=GridConfig(20000, 8000))
    for key in ("f90_hz", "f95_hz"):
        assert abs(fine[key] / coarse[key] - 1) < 0.02


def test_upper_band_edge_decreases_with_stellar_size():
    f90 = [_band(star=s)["f90_hz"] for s in (0.01, 0.02, 0.03, 0.04)]
    assert all(b < a for a, b in zip(f90, f90[1:])), f90


def test_quantiles_are_ordered():
    r = _band()
    assert r["f05_hz"] < r["f50_hz"] < r["f90_hz"] < r["f95_hz"] < r["f99_hz"]
    assert r["chord_time_s"] > 0


def test_interpolated_quantile_lies_between_bins():
    freq = np.arange(10, dtype=float)
    power = np.ones(10)
    q_int = cumulative_quantile_frequency(freq, power, 0.55, interpolate=True)
    q_bin = cumulative_quantile_frequency(freq, power, 0.55, interpolate=False)
    assert q_bin == 5.0 and 4.0 < q_int < 6.0 and q_int != q_bin


def test_baseline_one_keeps_dip_area_in_lowest_bins():
    x = np.linspace(-1000, 1000, 2001)
    dip = 1.0 - 0.5 * np.exp(-0.5 * (x / 20.0) ** 2)
    f, p1 = spatial_power_spectrum(x, dip, baseline=1.0)
    _, pm = spatial_power_spectrum(x, dip, baseline=None)
    c1, cm = np.cumsum(p1) / p1.sum(), np.cumsum(pm) / pm.sum()
    # mean subtraction zeroes the DC bin, i.e. erases the dip's real area term
    assert pm[0] < 1e-20 and p1[0] > 1e-3
    assert c1[2] > cm[2]
