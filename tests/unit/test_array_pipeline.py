"""
End-to-end tests of the 3-telescope pipeline on synthetic data: the
array injection Monte Carlo recovers an obvious event in coincidence,
and search_observation finds a simultaneous injected event in a
synthetic 3-channel observation while an event in a single channel is
rejected by the coincidence condition.
"""

import numpy as np
import pytest

from kbo_occultation import config
from kbo_occultation.config import (BandpassConfig, GridConfig, KBOConfig,
                                    NumericalConfig, ORM_ARRAY, StarConfig)
from kbo_occultation.simulation.detectability import spatial_to_time
from kbo_occultation.analysis.injection import inject_occultation
from kbo_occultation.analysis.injection_batch import efficiency_curve, run_array_injection_monte_carlo
from kbo_occultation.observations.photometry import LightCurve
from kbo_occultation.analysis.search import search_observation
from kbo_occultation.simulation.engine import compute_lightcurve

DT_S = config.standard_sampling * config.standard_sample_duration / 1e9
V_REL = 25_000.0
DISTANCE_AU = 43.0

STAR = StarConfig(temperature_K=6000.0, angular_radius_mas=0.0)  # point source (fast)
BAND = BandpassConfig(450.0, 550.0, 3)
GRID = GridConfig(x_max_m=5000.0, n_x=1000)
NUMERICS = NumericalConfig(n_int=40, n_r_grid=500, n_star_side=1)

DATA_FORMAT = [('magic_word', '<Q'),
               ('mean_chA', '<f'), ('mean_chB', '<f'), ('mean_chC', '<f'), ('mean_chD', '<f'),
               ('std_chA', '<f'), ('std_chB', '<f'), ('std_chC', '<f'), ('std_chD', '<f'),
               ('min_chA', '<f'), ('min_chB', '<f'), ('min_chC', '<f'), ('min_chD', '<f'),
               ('max_chA', '<f'), ('max_chB', '<f'), ('max_chC', '<f'), ('max_chD', '<f'),
               ('time_stamp', '<Q')]


def synthetic_channels(n=80_000, noise=0.05, seed=0):
    """Three uncorrelated flat-ish light curves on a shared time grid."""
    rng = np.random.default_rng(seed)
    t = DT_S * np.arange(n)
    lcs = {}
    for ch in ("A", "B", "C"):
        flux = 10.0 * (1.0 + noise * rng.standard_normal(n))
        lcs[ch] = LightCurve(t, flux, meta={"channel": ch})
    return lcs


@pytest.fixture(scope="module")
def template():
    kbo = KBOConfig(radius_m=1000.0, distance_au=DISTANCE_AU)
    return compute_lightcurve(kbo, STAR, BAND, GRID, NUMERICS)


def test_array_mc_recovers_obvious_event(template):
    lcs = synthetic_channels()
    identity = lambda lc: lc

    trials, profiles = run_array_injection_monte_carlo(
        lcs, ORM_ARRAY,
        radii_m=[1000.0],
        n_trials_per_radius=3,
        grid=GRID, numerics=NUMERICS, star=STAR, band=BAND, response=None,
        kinematics=V_REL, distance_au=DISTANCE_AU,
        method_specs={"raw": identity},
        pad_s=9.0,
        min_ntel=2,
        seed=1,
    )

    assert len(trials) == 3
    for tr in trials:
        assert tr.n_tel_recovered == 3, f"trial at t0={tr.t0_s} recovered {tr.recovered_per_channel}"
        assert tr.coincidence_recovered
        assert tr.combined_snr > tr.snr_per_channel["A"]

    # noise profiles: one per (method, channel, radius)
    assert len(profiles) == 3
    assert {p.channel for p in profiles} == {"A", "B", "C"}

    eff = efficiency_curve(trials, "raw")
    assert eff(1000.0) == 1.0
    assert eff(50.0) == 0.0  # outside measured range


def write_bin(path, flux_per_channel):
    n = len(next(iter(flux_per_channel.values())))
    rec = np.zeros(n, dtype=DATA_FORMAT)
    for ch, flux in flux_per_channel.items():
        rec[f"std_ch{ch}"] = np.sqrt(flux)  # loader squares std -> flux
    rec["time_stamp"] = (1e6 * (1000.0 + DT_S * np.arange(n))).astype(np.uint64)
    rec.tofile(str(path))
    return str(path)


def test_search_observation_finds_coincident_event(tmp_path, template):
    x_m, intensity = template
    n = 80_000
    rng = np.random.default_rng(5)
    t = DT_S * np.arange(n)
    t0 = t[0] + 1000.0 + 20.0  # loader reconstructs time from t0=1000s

    flux = {}
    for ch in ("A", "B", "C"):
        base = 10.0 * (1.0 + 0.05 * rng.standard_normal(n))
        inj = inject_occultation(t + 1000.0, base, x_m, intensity, V_REL, t0_s=t0)
        flux[ch] = inj.injected_signal
    fn = write_bin(tmp_path / "coincident.bin", flux)

    result = search_observation(
        fn, template_radii_m=(1000.0,), distance_au=DISTANCE_AU,
        star=STAR, band=BAND, grid=GRID, numerics=NUMERICS,
        kinematics=V_REL, min_ntel=3, n_background_shifts=20,
    )

    assert len(result.coincidences) >= 1
    ev = result.coincidences[0]
    assert ev.n_tel == 3
    assert abs(ev.time_s - t0) < result.tolerance_s
    assert result.accidental is not None


def test_search_observation_vetoes_single_telescope_event(tmp_path, template):
    x_m, intensity = template
    n = 80_000
    rng = np.random.default_rng(6)
    t = DT_S * np.arange(n)
    t0 = t[0] + 1000.0 + 20.0

    flux = {}
    for ch in ("A", "B", "C"):
        base = 10.0 * (1.0 + 0.05 * rng.standard_normal(n))
        if ch == "B":  # only LST-1 sees a "dip" (e.g. an instrumental artifact)
            inj = inject_occultation(t + 1000.0, base, x_m, intensity, V_REL, t0_s=t0)
            base = inj.injected_signal
        flux[ch] = base
    fn = write_bin(tmp_path / "single.bin", flux)

    result = search_observation(
        fn, template_radii_m=(1000.0,), distance_au=DISTANCE_AU,
        star=STAR, band=BAND, grid=GRID, numerics=NUMERICS,
        kinematics=V_REL, min_ntel=2, n_background_shifts=0,
    )

    # channel B fires alone -> no coincidence
    assert any(len(res.candidates) > 0 for res in result.per_telescope.values())
    assert result.coincidences == []
