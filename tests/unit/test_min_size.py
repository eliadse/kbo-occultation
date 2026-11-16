import numpy as np
import pytest

from kbo_occultation.analysis import min_size as ms
from kbo_occultation.analysis.matched_filter import sliding_matched_filter_snr

DT = 0.005
FS = 1.0 / DT
NPERSEG = 4096


def _grid():
    return np.linspace(0.0, FS / 2.0, NPERSEG // 2 + 1)


def _template(width_s=0.1, depth=0.05):
    t = np.arange(-0.5, 0.5, DT)
    return t, -depth * np.exp(-0.5 * (t / width_s) ** 2) * np.cos(8 * t) ** 2


def test_star_radius_from_mag():
    assert ms.star_radius_from_mag(7.0) == pytest.approx(0.03)
    assert ms.star_radius_from_mag(8.0) == pytest.approx(0.03 * 10 ** -0.2)
    assert ms.star_radius_from_mag(9.0) == pytest.approx(0.0119, rel=2e-2)
    targets = ms.generic_targets()
    assert targets[0].mag_B == 7.0 and targets[-1].mag_B == 9.0 and len(targets) == 9
    assert all(a.angular_radius_mas > b.angular_radius_mas for a, b in zip(targets, targets[1:]))


def test_catalogue_fills_missing_temperature_and_radius(tmp_path):
    path = tmp_path / "cat.csv"
    path.write_text("id,B,spt,rad\nstarA,7.5,B9V,\nstarB,8.5,,0.02\n")
    a, b = ms.load_target_catalogue(str(path), {"name": "id", "mag_B": "B", "spectral_type": "spt",
                                               "angular_radius_mas": "rad"})
    assert a.name == "starA" and a.angular_radius_mas == pytest.approx(ms.star_radius_from_mag(7.5))
    assert 9000 < a.temperature_K < 12000                      # from the B9 spectral type
    assert b.angular_radius_mas == 0.02 and b.temperature_K == 10000.0
    with pytest.raises(ValueError):
        ms.load_target_catalogue(str(path), {"mag_B": "nope"})


def test_expected_snr_white_noise_is_template_norm_over_sigma():
    _, tmpl = _template()
    sigma = 0.01
    psd = np.full(NPERSEG // 2 + 1, 2 * sigma ** 2 / FS)
    snr = ms.expected_snr(tmpl, _grid(), psd, FS)
    assert snr == pytest.approx(np.linalg.norm(tmpl) / sigma, rel=1e-6)


def test_expected_snr_matches_sliding_matched_filter_on_white_noise():
    tt, tmpl = _template(depth=0.1)
    sigma, n = 0.01, 20000
    psd = np.full(NPERSEG // 2 + 1, 2 * sigma ** 2 / FS)
    expected = ms.expected_snr(tmpl, _grid(), psd, FS)
    time = np.arange(n) * DT
    peaks = []
    for seed in range(8):
        rng = np.random.default_rng(seed)
        signal = 1.0 + sigma * rng.standard_normal(n)
        i0 = 10000
        signal[i0:i0 + len(tmpl)] *= 1.0 + tmpl
        res = sliding_matched_filter_snr(time, signal, tt, 1.0 + tmpl, normalization="global")
        peaks.append(res.snr[np.argmin(np.abs(res.time - time[i0 + len(tmpl) // 2]))])
    assert np.mean(peaks) == pytest.approx(expected, rel=0.08)


def test_corrections_ordering_on_red_noise():
    _, tmpl = _template(width_s=0.2)
    f = _grid()
    psd = 2 * 0.002 ** 2 / FS + 1e-4 / np.maximum(f, 0.05) ** 2
    snr = {name: ms.expected_snr(tmpl, f, psd, FS, ms.correction_gain(name, f, psd, FS))
           for name in ("raw", "whitened", "optimal", "butter_1")}
    assert snr["optimal"] >= snr["whitened"] > snr["raw"]
    assert snr["optimal"] >= snr["butter_1"] > snr["raw"]


def test_butterworth_gain_is_zero_phase_squared_response():
    f = _grid()
    g = ms.correction_gain("butter_1", f, np.ones_like(f), FS)
    assert g[0] == pytest.approx(0.0, abs=1e-12)
    assert g[np.searchsorted(f, 1.0)] == pytest.approx(0.5, rel=0.12)     # forward-backward: |H|^2, = 1/2 at fc
    assert g[-1] == pytest.approx(1.0, abs=1e-3)


def test_detection_probability_matches_monte_carlo():
    rng = np.random.default_rng(0)
    snr = {"A": np.array([2.0, 5.0, 8.0]), "B": np.array([3.0, 6.0, 9.0]), "C": np.array([1.0, 4.0, 7.0])}
    thr, scatter = 5.0, 1.3
    draws = {k: v[None, :] + scatter * rng.standard_normal((200000, 3)) for k, v in snr.items()}
    n_above = sum((d > thr).astype(int) for d in draws.values())
    for k in (1, 2, 3):
        p = ms.detection_probability(snr, thr, k, scatter)
        assert p == pytest.approx((n_above >= k).mean(axis=0), abs=0.01)


def test_r50_interpolates_in_log_radius_and_flags():
    radii = np.array([10.0, 100.0, 1000.0])
    assert ms.r50(radii, [0.0, 0.5, 1.0]) == (pytest.approx(100.0), "ok")
    value, status = ms.r50(radii, [0.0, 0.25, 0.75])
    assert status == "ok" and value == pytest.approx(100.0 * 10 ** 0.5)
    assert ms.r50(radii, [0.7, 0.8, 0.9]) == (10.0, "below_grid")
    assert np.isnan(ms.r50(radii, [0.0, 0.1, 0.2])[0])
    # a dip back below 0.5 at a larger radius moves R50 up (conservative)
    assert ms.r50(np.array([1.0, 2.0, 3.0, 4.0]), [0.6, 0.4, 0.6, 0.9])[0] > 2.0


def test_rate_curve_empirical_and_exponential_tail():
    rng = np.random.default_rng(1)
    hours, lam, scale = 100.0, 40.0, 1.0
    n = int(lam * hours)
    snr = 3.0 + rng.exponential(scale, n)
    curve = ms.RateCurve(snr, hours, n_tail=200)
    assert curve.rate(3.5) == pytest.approx(lam * np.exp(-0.5), rel=0.05)
    assert not curve.extrapolated(3.5)
    s = 3.0 + 8.0
    assert curve.extrapolated(s)
    assert curve.rate(s) == pytest.approx(lam * np.exp(-8.0), rel=0.4)
    assert ms.RateCurve(np.array([]), 1.0).rate(5.0) == 0.0


def test_coincidence_threshold_matches_closed_form():
    rng = np.random.default_rng(2)
    hours, lam, scale = 200.0, 40.0, 1.0
    curves = {ch: ms.RateCurve(3.0 + rng.exponential(scale, int(lam * hours)), hours, n_tail=1000)
              for ch in "ABC"}
    tol, far = 0.4, 0.001
    s, extrapolated = ms.coincidence_threshold(curves, tol, far, min_ntel=2)
    r = np.sqrt(far / (3 * 2 * tol / 3600.0))                  # 3 pairs, 2 tau r^2 each
    assert s == pytest.approx(3.0 - scale * np.log(r / lam), abs=0.15)
    assert extrapolated
    s3, _ = ms.coincidence_threshold(curves, tol, far, min_ntel=3)
    assert s3 < s                                              # the 3-fold coincidence is rarer
    s1, extrap1 = ms.single_threshold(curves["A"], 1.0)
    assert s1 == pytest.approx(3.0 + scale * np.log(lam), abs=0.15) and extrap1
    with pytest.raises(ValueError):
        ms.coincidence_threshold({"A": curves["A"]}, tol, far, min_ntel=2)


def test_coincidence_rate_is_quadratic_in_the_single_rate():
    flat = {ch: ms.RateCurve(np.array([]), 1.0) for ch in "AB"}
    assert ms.coincidence_rate(flat, 5.0, 0.4) == 0.0


def test_profiles_roundtrip(tmp_path):
    f = _grid()
    p = ms.NoiseProfile("n", "o", "A", "MAGIC-2", 7.0, None, DT, f, np.abs(np.sin(f)), 0.9, 0.95, "measured",
                        1.2, 0.9, 1.0, 0.23, 0.1, 30.0)
    path = str(tmp_path / "p.npz")
    ms.save_profiles([p], path)
    q, = ms.load_profiles(path)
    assert q.telescope == "MAGIC-2" and q.teff_K is None and q.fs == pytest.approx(FS)
    np.testing.assert_allclose(q.S_excess, p.S_excess)


def _profile(telescope, channel):
    f = _grid()
    return ms.NoiseProfile("n", "o", channel, telescope, 7.5, None, 10 * ms.FRAME_S, f,
                           np.zeros_like(f), 0.9, 0.9, "model", 1.2, 1.0, 1.0, 0.23, 0.0, 30.0)


@pytest.fixture(scope="module")
def tiny_cfg():
    return ms.MinSizeConfig(radii_m=np.geomspace(30.0, 600.0, 7), x_half_m=6000.0, dx_m=8.0,
                            n_star_side=16, n_r_grid=1500)


def test_theory_r50_grows_with_magnitude_and_shrinks_with_lst_focus(tiny_cfg):
    cfg = tiny_cfg
    profile = _profile("LST-1", "B")
    curves = {"B": ms.RateCurve(3.0 + np.random.default_rng(3).exponential(1.0, 800), 20.0)}
    targets = [ms.TargetStar(f"B{m}", m, 0.03) for m in (9.0, 11.0, 13.0)]    # fixed size: only the noise changes
    table = ms.min_size_table(targets, [profile], cfg, curves, noise_models=("theory",), focus_gains=(1.0, 2.0))
    t = table[(table.config == "LST-1") & (table.far_per_hr == cfg.coinc_far_per_hr)]
    r_plain = t[t.lst_focus == 1.0].sort_values("mag_B")["r50_m"].to_numpy()
    r_focus = t[t.lst_focus == 2.0].sort_values("mag_B")["r50_m"].to_numpy()
    assert np.all(np.isfinite(r_plain)) and np.all(np.diff(r_plain) > 0)
    assert np.all(r_focus <= r_plain)
    assert set(table.far_per_hr) == {1.0, 0.1, 0.01}
    assert (ms.summarize(table)["n_runs"] == 1).all()


def test_excess_noise_raises_r50(tiny_cfg):
    cfg = tiny_cfg
    f = _grid()
    clean, noisy = _profile("MAGIC-2", "A"), _profile("MAGIC-2", "A")
    noisy.S_excess = 1e-5 / np.maximum(f, 0.1) ** 2
    curves = {"A": ms.RateCurve(3.0 + np.random.default_rng(4).exponential(1.0, 800), 20.0)}
    star = ms.generic_targets([8.0])[0]
    r = {}
    for name, prof in (("clean", clean), ("noisy", noisy)):
        table = ms.min_size_table([star], [prof], cfg, curves, noise_models=("measured",))
        r[name] = table[(table.config == "MAGIC-2") & (table.far_per_hr == 1.0)]["r50_m"].iloc[0]
    assert r["noisy"] > r["clean"]


# ───────────────────────────── calibration against injections ───────────────

import pandas as pd


def _thresholds(single=6.0, coinc=4.0):
    rows = []
    for far in (1.0,):
        rows += [{"far_per_hr": far, "config": ch, "k": 1, "threshold": single} for ch in "ABC"]
        rows += [{"far_per_hr": far, "config": "2of3", "k": 2, "threshold": coinc},
                 {"far_per_hr": far, "config": "3of3", "k": 3, "threshold": coinc}]
    return pd.DataFrame(rows)


def _hand_made(snr_by_channel, t_by_channel, t0=100.0):
    """One run, one radius, one trial with the given candidate SNRs / times per channel."""
    key = {"night": "n", "observation": "o", "correction": "raw", "sf_mode": "measured"}
    trials = pd.DataFrame([{**key, "radius_m": 100.0, "trial": 0, "channel": ch, "t0_s": t0, "peak_snr": 9.0}
                           for ch in "ABC"])
    cands = pd.DataFrame([{**key, "radius_m": 100.0, "trial": 0, "channel": ch, "t0_s": t0, "time_s": t_by_channel[ch],
                           "dt_s": t_by_channel[ch] - t0, "snr": snr, "chi2_reduced": 1.0, "template_radius_m": 100.0}
                          for ch, snr in snr_by_channel.items()])
    return trials, cands


def _eff(eff, config):
    return float(eff[(eff["config"] == config)]["eff"].iloc[0])


def test_measured_efficiency_k_of_n_and_thresholds():
    # A and B above the coincidence threshold (4) and close in time, C absent
    trials, cands = _hand_made({"A": 5.0, "B": 4.5}, {"A": 100.0, "B": 100.1})
    eff = ms.measured_efficiency(trials, cands, _thresholds(), tolerance_s=0.4)
    assert _eff(eff, "2of3") == 1.0 and _eff(eff, "3of3") == 0.0
    assert _eff(eff, "A") == 0.0 and _eff(eff, "B") == 0.0          # below the single threshold (6)
    # below the coincidence threshold: ignored
    trials, cands = _hand_made({"A": 5.0, "B": 3.9}, {"A": 100.0, "B": 100.1})
    assert _eff(ms.measured_efficiency(trials, cands, _thresholds(), 0.4), "2of3") == 0.0
    # too far apart in time: no coincidence event
    trials, cands = _hand_made({"A": 5.0, "B": 5.0}, {"A": 100.0, "B": 100.6})
    assert _eff(ms.measured_efficiency(trials, cands, _thresholds(), 0.4), "2of3") == 0.0
    # single channel above its threshold
    trials, cands = _hand_made({"A": 7.0}, {"A": 100.0})
    eff = ms.measured_efficiency(trials, cands, _thresholds(), 0.4)
    assert _eff(eff, "A") == 1.0 and _eff(eff, "B") == 0.0 and _eff(eff, "2of3") == 0.0


def test_excess_scaling_slopes():
    rng = np.random.default_rng(0)
    f = _grid()
    for slope, expect in ((0.0, "star"), (-2.0, "total")):
        profiles = []
        for tel, ch in (("MAGIC-1", "C"), ("MAGIC-2", "A"), ("LST-1", "B")):
            for sf in rng.uniform(0.3, 0.95, 8):
                p = _profile(tel, ch)
                p.star_fraction = float(sf)
                p.S_excess = 1e-4 * sf ** slope / np.maximum(f, 0.5) ** 2 * (1.0 + 0.02 * rng.standard_normal(len(f)))
                profiles.append(p)
        out = ms.excess_scaling_test(profiles).set_index("group")
        assert out.loc["pooled", "slope"] == pytest.approx(slope, abs=0.1), expect
        assert out.loc["LST", "n"] == 8


def test_fit_snr_calibration_recovers_scale():
    rng = np.random.default_rng(1)
    key = {"night": "n", "observation": "o", "correction": "raw", "sf_mode": "measured"}
    radii = [50.0, 100.0, 200.0]
    pred = pd.DataFrame([{"night": "n", "observation": "o", "correction": "raw", "channel": "A",
                          "telescope": "MAGIC-2", "radius_m": r, "snr_pred": s} for r, s in zip(radii, (4.0, 8.0, 16.0))])
    trials = pd.DataFrame([{**key, "channel": "A", "radius_m": r, "trial": i, "peak_snr": 1.5 * s + rng.standard_normal()}
                           for r, s in zip(radii, (4.0, 8.0, 16.0)) for i in range(40)])
    per_run, summary = ms.fit_snr_calibration(trials, pred)
    assert per_run["scale"].iloc[0] == pytest.approx(1.5, rel=0.05)
    assert per_run["scatter"].iloc[0] == pytest.approx(1.0, rel=0.4)
    assert summary[(summary["kind"] == "all")]["scale"].iloc[0] == pytest.approx(1.5, rel=0.05)


def test_predicted_efficiency_scale_shifts_the_curve():
    pred = pd.DataFrame([{"night": "n", "observation": "o", "correction": "raw", "channel": ch, "telescope": "x",
                          "radius_m": r, "snr_pred": s} for ch in "ABC" for r, s in zip((50.0, 100.0), (4.0, 8.0))])
    e1 = ms.predicted_efficiency(pred, _thresholds())
    e2 = ms.predicted_efficiency(pred, _thresholds(), {"raw": (1.5, 1.0)})
    a = e1[e1.config == "2of3"].sort_values("radius_m")["eff_pred"].to_numpy()
    b = e2[e2.config == "2of3"].sort_values("radius_m")["eff_pred"].to_numpy()
    assert np.all(b >= a) and b[0] > a[0] + 0.1


def test_scale_by_deviance_finds_true_scale():
    pred = pd.DataFrame([{"night": "n", "observation": "o", "correction": "raw", "channel": ch, "telescope": "x",
                          "radius_m": r, "snr_pred": 1.2 * r / 10.0} for ch in "ABC" for r in (20.0, 25.0, 30.0, 35.0, 40.0)])
    thr = _thresholds(single=6.0, coinc=6.0)
    truth = ms.predicted_efficiency(pred, thr, {"raw": (1.3, 1.0)})
    meas = truth[truth["config"] == "2of3"].copy()
    meas["n"] = 400
    meas["hits"] = (meas["eff_pred"] * 400).round().astype(int)
    meas["sf_mode"] = "measured"
    meas = meas.drop(columns="eff_pred")
    best, curve = ms.scale_by_deviance(meas, pred, thr, "raw")
    assert best == pytest.approx(1.3, rel=0.08)


def test_timeslide_far_matches_independent_streams(tmp_path):      # sparse lists: window occupancy ~ 5 %
    from types import SimpleNamespace
    rng = np.random.default_rng(5)
    t0 = 1.7e9
    live = 300.0
    rows = []
    for ch in "ABC":
        for t, s in zip(t0 + 10 + rng.uniform(0, live - 20, 40), 3.0 + rng.exponential(0.8, 40)):
            rows.append({"night": "n", "observation": "o", "channel": ch, "method": "raw", "time_s": t, "snr": s,
                         "chi2_reduced": 1.0, "radius_m": 100.0, "depth": 0.01})
    csv = str(tmp_path / "c.csv")
    pd.DataFrame(rows).to_csv(csv, index=False)
    tel = lambda ch: SimpleNamespace(name=ch)                                                  # noqa: E731
    obs = SimpleNamespace(name="o", has_cache=True, good_time=[0.0, live], time_span=lambda: (t0, t0 + live))
    night = SimpleNamespace(name="n", channels=["A", "B", "C"], stars={"o": obs}, telescope=tel)
    log = SimpleNamespace(select=lambda nights=None, include_excluded=False: [night])
    curves = {ch: ms.RateCurve(np.array([r["snr"] for r in rows if r["channel"] == ch]), live / 3600.0)
              for ch in "ABC"}
    tol = 0.2
    far = ms.timeslide_far(csv, log, "raw", tolerance_s=tol, snr_grid=[3.0, 3.5], n_shifts=400, curves=curves)
    for _, r in far.iterrows():
        assert r["slide_rate_per_hr"] == pytest.approx(r["analytic_rate_per_hr"], rel=0.25)
    assert far["slide_hours"].iloc[0] == pytest.approx(400 * live / 3600.0)


def test_calibrated_config_and_theory_noise_is_not_scaled(tiny_cfg):
    from dataclasses import replace
    c = replace(tiny_cfg, correction="whitened").calibrated()
    assert (c.snr_scale, c.snr_scatter) == ms.SNR_CALIBRATION["whitened"]
    unknown = replace(tiny_cfg, correction="nonexistent", snr_scale=1.3).calibrated()
    assert unknown.snr_scale == 1.3
    star, profile = ms.TargetStar("t", 9.0, 0.03), _profile("MAGIC-2", "A")
    templates = ms.template_family(star, tiny_cfg, __import__("kbo_occultation.simulation.instruments",
                                                                 fromlist=["Instrument"]).Instrument("MAGIC", site="LaPalma"))
    base = ms.per_telescope_snr(star, profile, templates, tiny_cfg, theory_only=True)
    scaled = ms.per_telescope_snr(star, profile, templates, replace(tiny_cfg, snr_scale=2.0), theory_only=True)
    np.testing.assert_allclose(base, scaled)
    assert np.allclose(2.0 * ms.per_telescope_snr(star, profile, templates, tiny_cfg),
                       ms.per_telescope_snr(star, profile, templates, replace(tiny_cfg, snr_scale=2.0)))
