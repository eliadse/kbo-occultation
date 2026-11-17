import numpy as np
import pandas as pd
import pytest
from scipy.signal import welch

from kbo_occultation.analysis import min_size as ms
from kbo_occultation.analysis import synthetic as sy
from kbo_occultation.analysis.noise_removal import ComparisonConfig

FS = 1.0 / ms.FRAME_S


@pytest.fixture(scope="module")
def cfg():
    return ms.MinSizeConfig(time_resolved=False)


def test_gaussian_noise_white_level_survives_binning():
    sigma_bin = 0.01
    fs_bin = FS / 10
    x = sy.gaussian_noise(1_000_000, FS, 2.0 * sigma_bin ** 2 / fs_bin, rng=np.random.default_rng(0))
    assert x.shape == (1_000_000,)
    binned = x.reshape(-1, 10).mean(axis=1)
    assert binned.std() == pytest.approx(sigma_bin, rel=0.02)
    assert x.std() == pytest.approx(sigma_bin * np.sqrt(10), rel=0.02)


def test_gaussian_noise_coloured_matches_the_requested_psd():
    fs_bin = FS / 10
    white = 2.0 * 0.01 ** 2 / fs_bin
    freq = np.linspace(0.0, fs_bin / 2.0, 2049)
    excess = 1e-4 / np.maximum(freq, 0.1) ** 2
    x = sy.gaussian_noise(2_000_000, FS, white, freq, excess, np.random.default_rng(1))
    binned = x.reshape(-1, 10).mean(axis=1)
    f, psd = welch(binned, fs=fs_bin, nperseg=4096)
    for lo, hi in ((0.2, 0.5), (1.0, 3.0), (5.0, 20.0), (50.0, 90.0)):
        sel = (f >= lo) & (f <= hi)
        target = np.interp(f[sel], freq, excess) + white
        assert np.mean(psd[sel] / target) == pytest.approx(1.0, rel=0.1)


def test_synthetic_profile_validates_and_sets_the_white_sigma(cfg):
    with pytest.raises(ValueError):
        sy.synthetic_profile("n", "o", "A", 8.0, cfg, noise="coloured")           # no reference
    with pytest.raises(ValueError):
        sy.synthetic_profile("n", "o", "A", 8.0, cfg, noise="pink")
    with pytest.raises(ValueError):
        sy.synthetic_profile("n", "o", "A", 8.0, cfg, snr=10.0)                    # no snr_dt_s
    p = sy.synthetic_profile("n", "o", "B", 8.0, cfg, snr=10.0, snr_dt_s=4.0 * cfg.bin_dt_s)
    assert p.telescope == "LST-1" and p.star_fraction == 1.0
    assert p.white_sigma == pytest.approx(0.1 * np.sqrt(4.0))                      # 10 per 4 bins -> sigma per bin
    assert np.all(p.S_excess == 0.0)


def test_noise_psd_honours_white_sigma_and_it_survives_save_load(cfg, tmp_path):
    star = ms.TargetStar("s", 8.0, 0.02)
    p = sy.synthetic_profile("n", "o", "A", 8.0, cfg, snr=10.0, snr_dt_s=cfg.bin_dt_s)
    assert ms.noise_psd(p, star, cfg)[0] == pytest.approx(2.0 * 0.1 ** 2 / p.fs)
    model = sy.synthetic_profile("n", "o", "A", 8.0, cfg)
    assert model.white_sigma is None
    path = str(tmp_path / "p.npz")
    ms.save_profiles([p, model], path)
    a, b = ms.load_profiles(path)
    assert a.white_sigma == pytest.approx(p.white_sigma) and b.white_sigma is None


def test_observed_target_uses_the_synthetic_radius_and_keeps_the_log_default():
    class Logged:
        name, mag_B, spectral_type = "HD1", 7.5, None

    class Synthetic(Logged):
        angular_radius_mas = 0.0189

    assert ms.observed_target(Logged()).angular_radius_mas == 0.03
    assert ms.observed_target(Synthetic()).angular_radius_mas == pytest.approx(0.0189)


def test_synthetic_run_is_deterministic_independent_between_channels_and_has_the_model_level(cfg):
    night, obs, group = sy.synthetic_runs([8.0], cfg, n_runs=1, run_s=120.0, seed=3)[0]
    assert night.name == "synthetic-white" and obs.name == "B8_0" and [p.channel for p in group] == ["A", "B", "C"]
    assert obs.angular_radius_mas == pytest.approx(ms.star_radius_from_mag(8.0))
    t0, t1 = obs.time_span()
    assert t1 - t0 == pytest.approx(120.0, abs=0.01)
    a1, a2, b = obs.lightcurve("A"), obs.lightcurve("A"), obs.lightcurve("B")
    np.testing.assert_array_equal(a1.signal, a2.signal)
    a2.signal[:] = 0.0                                                       # a copy: the cache is not touched
    np.testing.assert_array_equal(obs.lightcurve("A").signal, a1.signal)
    assert abs(np.corrcoef(a1.signal, b.signal)[0, 1]) < 0.01
    sigma_bin = np.sqrt(ms.noise_psd(group[0], ms.observed_target(obs), cfg)[0] * group[0].fs / 2.0)
    assert a1.signal.mean() == pytest.approx(1.0, abs=1e-3)
    assert a1.signal.std() == pytest.approx(sigma_bin * np.sqrt(cfg.bin_samples), rel=0.03)
    obs.release()
    np.testing.assert_array_equal(obs.lightcurve("A").signal, a1.signal)     # redrawn identically


def test_reference_runs_pick_the_median_excess_and_coloured_runs_carry_it(cfg):
    def profile(night, obs, ch, level):
        p = sy.synthetic_profile(night, obs, ch, 8.0, cfg)
        p.S_excess = level * 1e-4 / np.maximum(p.freq, 0.1) ** 2
        return p

    profiles = [profile("n", f"star{i}", ch, lvl) for i, lvl in enumerate((1.0, 5.0, 2.0)) for ch in "ABC"]
    profiles.append(profile("n", "partial", "A", 100.0))                      # not a complete run
    assert [g["A"].observation for g in sy.reference_runs(profiles, "median")] == ["star2"]
    assert [g["A"].observation for g in sy.reference_runs(profiles, "all")] == ["star0", "star2", "star1"]
    assert [g["A"].observation for g in sy.reference_runs(profiles, "n/star1")] == ["star1"]
    with pytest.raises(ValueError):
        sy.reference_runs(profiles, "n/partial")
    runs = sy.synthetic_runs([8.0], cfg, noise="coloured", references=sy.reference_runs(profiles, "n/star1"),
                             n_runs=1, run_s=60.0)
    assert runs[0][0].name == "synthetic-coloured"
    np.testing.assert_allclose(runs[0][2][0].S_excess, profiles[3].S_excess)
    with pytest.raises(ValueError):
        sy.synthetic_runs([8.0], cfg, noise="coloured", n_runs=1)


def test_rate_curves_from_candidates_pool_hours_and_counts():
    calib = pd.DataFrame([{"night": "n", "observation": f"o{i}", "correction": "raw", "sf_mode": sf, "channel": ch,
                           "live_s": 1800.0} for i in range(2) for sf in ("measured", "model") for ch in "AB"])
    cands = pd.DataFrame([{"night": "n", "observation": "o0", "correction": "raw", "channel": "A", "snr": s}
                          for s in (3.5, 4.5, 6.0)])
    curves = sy.rate_curves_from_candidates(cands, calib)
    assert set(curves) == {("n", "raw")} and set(curves[("n", "raw")]) == {"A", "B"}
    a, b = curves[("n", "raw")]["A"], curves[("n", "raw")]["B"]
    assert a.hours == pytest.approx(1.0) and len(a.snr) == 3 and len(b.snr) == 0       # model-mode rows not doubled
    assert a.rate(3.5) == pytest.approx(2.0) and b.rate(3.0) == 0.0          # empirical count up to the lowest of 3


def test_end_to_end_peak_snr_matches_the_search_bank_prediction(cfg):
    ccfg = ComparisonConfig(n_trials=3, seed=0)
    radii = [100.0, 200.0]
    runs = sy.synthetic_runs([8.0], cfg, n_runs=1, run_s=150.0, seed=1)
    night, obs, group = runs[0]
    full = []
    trials, cands, calib = ms.calibrate_run(night, obs, group, cfg, ccfg, "whitened", radii, n_trials=3,
                                            full_candidates=full)
    assert len(trials) == 3 * 3 * len(radii) and set(trials["channel"]) == {"A", "B", "C"}
    assert len(full) > 0 and {"channel", "telescope", "snr", "night", "observation"} <= set(full[0])
    assert np.all(np.abs(calib["snr_std"] - 1.0) < 0.15)                      # empirical normalisation on white noise
    from kbo_occultation.analysis.mimic import observed_star_mas

    own_bank = ms.bank_residuals(cfg, star=ms.observed_target(obs, observed_star_mas(obs)))   # calibrate_run's default bank
    pred = ms.predicted_snr(runs, cfg, ["whitened"], radii, bank=own_bank)
    got = trials.groupby(["channel", "radius_m"])["peak_snr"].median()
    for row in pred.itertuples():
        assert got[(row.channel, row.radius_m)] == pytest.approx(row.snr_pred, rel=0.2)
    matched = ms.predicted_snr(runs, cfg, ["whitened"], radii)
    assert np.all(matched["snr_pred"].to_numpy() >= pred["snr_pred"].to_numpy() - 1e-9)   # the bank cannot beat a perfect match
