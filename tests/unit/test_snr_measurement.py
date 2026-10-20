import os

import astropy.units as u
import numpy as np
import pytest
from scipy.special import gammaln

from kbo_occultation.simulation.instruments import Instrument
from kbo_occultation.observations.log import load_observation_log
from kbo_occultation.observations.photometry import LightCurve
from kbo_occultation.analysis.snr_measurement import (detrended_snr, extinction_factor, good_time_lightcurve,
                                             histogram_peak, measured_snr, throughput_from_faulty_panels)

DT = 5.24288e-4


def white_lc(mean, sigma, n=200_000, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(n) * DT
    return LightCurve(t, rng.normal(mean, sigma, n))


def test_measured_and_detrended_snr_white_noise():
    # White noise: both estimators should give (mean - nsb) / (sigma / sqrt(n_bin)).
    lc = white_lc(50.0, 2.0)
    n_bin, nsb = 10, 8.0
    expected = (50.0 - nsb) / (2.0 / np.sqrt(n_bin))
    m = measured_snr(lc, n_bin, nsb)
    d = detrended_snr(lc, n_bin, nsb)
    # measured_snr averages per-bin std's (ddof=0): E[s] = sigma sqrt((n-1)/n) c4(n)
    c4 = np.sqrt(2 / (n_bin - 1)) * np.exp(gammaln(n_bin / 2) - gammaln((n_bin - 1) / 2))
    bias = np.sqrt((n_bin - 1) / n_bin) * c4
    assert m["snr"] == pytest.approx(expected / bias, rel=0.02)
    assert d["snr"] == pytest.approx(expected, rel=0.05)


def test_detrended_snr_ignores_slow_drift():
    lc = white_lc(50.0, 2.0)
    drifting = LightCurve(lc.time, lc.signal * (1 + 0.05 * np.sin(2 * np.pi * lc.time / 60.0)))
    flat = detrended_snr(lc, 10, 0.0)["snr"]
    assert detrended_snr(drifting, 10, 0.0)["snr"] == pytest.approx(flat, rel=0.05)


def test_histogram_peak_last_and_first():
    rng = np.random.default_rng(1)
    data = np.concatenate([rng.normal(10, 0.5, 5000), rng.normal(30, 0.5, 20000)])
    assert histogram_peak(data, "last") == pytest.approx(30, abs=0.5)
    assert histogram_peak(data, "first") == pytest.approx(10, abs=0.5)
    # mean peak is dragged down by the low population; peak-based signal is not
    lc = LightCurve(np.arange(data.size) * DT, data)
    assert measured_snr(lc, 10, 0.0, signal="peak")["signal"] > measured_snr(lc, 10, 0.0)["signal"]


def test_good_time_and_extinction_and_throughput():
    lc = white_lc(1.0, 0.1, n=10_000)
    sel = good_time_lightcurve(lc, [1.0, 2.0])
    assert sel.time[0] - lc.time[0] >= 1.0 and sel.time[-1] - lc.time[0] <= 2.0
    assert good_time_lightcurve(lc, None) is lc
    assert extinction_factor(0.0, 0.23) == pytest.approx(1.0)
    assert extinction_factor(60.0, 0.23) == pytest.approx(10 ** (0.4 * 0.23))
    assert throughput_from_faulty_panels(None, "m1") == 1.0
    assert throughput_from_faulty_panels(2, "m2") == pytest.approx(1 - 2 / 247)


def test_noise_budget_new_options_default_to_old_result():
    tel = Instrument("MAGIC", site="LaPalma")
    base = tel.noise_budget(7.0, 5 * u.ms, airmass=1.5)
    same = tel.noise_budget(7.0, 5 * u.ms, airmass=1.5, extinction=False, throughput_factor=1.0)
    assert same["sigma_total"] == pytest.approx(base["sigma_total"])
    # extinction at X=1.5 and lost mirror area both reduce the star photons
    ext = tel.noise_budget(7.0, 5 * u.ms, airmass=1.5, extinction=True)
    assert ext["star_photons"] < base["star_photons"]
    assert tel.airmass == 1.0  # restored after the call
    half = tel.noise_budget(7.0, 5 * u.ms, throughput_factor=0.5)
    assert half["star_photons"] == pytest.approx(0.5 * tel.noise_budget(7.0, 5 * u.ms)["star_photons"])


def test_observation_log_files_exist():
    log = load_observation_log()
    for obs in log.all_observations():
        assert os.path.exists(obs.path), obs.path
    assert log.nights["2026-03-08"].telescopes == ["MAGIC-2", "LST-1", "MAGIC-1"]


def test_frame_noise_white_and_estimator_floor():
    from kbo_occultation.analysis.snr_measurement import frame_noise, variance_estimator_floor
    lc = white_lc(50.0, 0.5)
    assert frame_noise(lc) == pytest.approx(0.01, rel=0.02)
    assert variance_estimator_floor(2**18) == pytest.approx(np.sqrt(2 / 2**18))
    # the variance of N Gaussian samples really has relative error sqrt(2/N)
    rng = np.random.default_rng(3)
    v = rng.normal(0, 1, (2000, 4096)).var(axis=1)
    assert np.std(v) == pytest.approx(variance_estimator_floor(4096), rel=0.05)
