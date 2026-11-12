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
