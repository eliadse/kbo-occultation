import pandas as pd
import numpy as np
import pytest
import astropy.units as u

from kbo_occultation.analysis.injection import inject_occultation
from kbo_occultation.analysis.noise_removal import (make_whitener, threshold_for_far, Method, _vetoed)
from kbo_occultation.config import DC_CUT_HZ
from kbo_occultation.observations.dc_combine import (check_dc_crossover, combine_fast_and_dc,
                                                     raised_cosine_lowpass)
from kbo_occultation.observations.photometry import LightCurve
from kbo_occultation.simulation.instruments import Instrument


def test_raised_cosine_is_half_at_crossover_and_monotonic():
    f = np.linspace(0, 1, 1001)
    w = raised_cosine_lowpass(f, 0.3, 0.5)
    assert w[0] == 1.0 and w[-1] == 0.0
    assert w[np.argmin(np.abs(f - 0.3))] == pytest.approx(0.5, abs=1e-2)
    assert np.all(np.diff(w) <= 1e-12)


def test_crossover_validation():
    check_dc_crossover(0.3, 0.5)
    with pytest.raises(ValueError):
        check_dc_crossover(1.0, 0.5)          # the old default, above the DC Nyquist
    with pytest.raises(ValueError):
        check_dc_crossover(0.4, 0.5, 0.5)     # transition would reach past Nyquist
    for key, cut in DC_CUT_HZ.items():
        check_dc_crossover(cut, 0.5 if key != "lst1" else 0.25)


def test_dc_swap_keeps_fast_above_crossover_and_dc_below():
    dt, n = 0.01, 2 ** 15
    t = np.arange(n) * dt
    rng = np.random.default_rng(0)
    fast = 1.0 + 0.01 * rng.standard_normal(n)
    dc = 1.0 + 0.05 * np.sin(2 * np.pi * 0.05 * t)       # slow trend only the DC knows about
    out = combine_fast_and_dc(t, fast, dc, 0.3, taper_s=5.0)
    mid = slice(n // 4, 3 * n // 4)
    # slow trend imported
    assert np.corrcoef(out[mid], dc[mid])[0, 1] > 0.95
    # ends fall back to the fast signal
    assert np.allclose(out[:20], fast[:20], atol=0.03)
    # high-frequency content untouched: differences are smooth
    resid = out - fast
    spec = np.abs(np.fft.rfft(resid))
    f = np.fft.rfftfreq(n, dt)
    assert spec[f > 0.6].max() < 0.02 * spec.max()


def test_dc_swap_identity_when_dc_equals_fast():
    t = np.arange(4096) * 0.01
    x = 2.0 + np.sin(t)
    assert np.allclose(combine_fast_and_dc(t, x, x, 0.3), x)


def test_star_fraction_leaves_background_undimmed():
    t = np.arange(0, 10, 0.001)
    x = np.array([-100.0, 0.0, 100.0])
    inten = np.array([1.0, 0.0, 1.0])
    full = inject_occultation(t, np.ones_like(t), x, inten, 1000.0, t0_s=5.0)
    half = inject_occultation(t, np.ones_like(t), x, inten, 1000.0, t0_s=5.0, star_fraction=0.5)
    assert full.injected_signal.min() == pytest.approx(0.0, abs=0.02)
    assert half.injected_signal.min() == pytest.approx(0.5, abs=0.02)


def test_noise_budget_new_terms():
    tel = Instrument("MAGIC")
    base = tel.noise_budget(6.0, 524.288 * u.us, variance_floor=False, scint_coeff=1.5)
    new = tel.noise_budget(6.0, 524.288 * u.us)
    # radiometer k = 250/110, scaled by the whole variance relative to the star: (S + B) / S
    ratio = (new["star_photons"] + new["nsb_photons"]) / new["star_photons"]
    assert new["sigma_floor"] == pytest.approx(np.sqrt(2 * (250 / 110) / 2 ** 18) * ratio, rel=1e-6)
    assert new["sigma_total"] > base["sigma_total"]
    assert base["sigma_floor"] == 0.0
    # the floor averages down with the number of frames in a bin, red noise does not
    coarse = tel.noise_budget(6.0, 5 * u.ms, red_noise_rms=2e-3)
    assert coarse["sigma_floor"] < new["sigma_floor"] / 3
    assert coarse["sigma_red"] == pytest.approx(2e-3)
    assert Instrument("LST").noise_budget(6.0, 524.288 * u.us)["sigma_floor"] > 0   # radiometer k applies to LST-1 too
    # default scintillation coefficient is 1.5
    s1 = tel.noise_budget(6.0, 5 * u.ms, scint_coeff=1.0)["sigma_scint"]
    assert tel.noise_budget(6.0, 5 * u.ms)["sigma_scint"] == pytest.approx(1.5 * s1)


def test_whitener_flattens_red_noise_and_keeps_mean():
    rng = np.random.default_rng(1)
    n, dt = 2 ** 18, 0.001
    white = rng.standard_normal(n)
    red = np.cumsum(rng.standard_normal(n)); red = (red - red.mean()) / red.std() * 3
    sig = 100.0 * (1 + 0.01 * (white + 0.5 * red))
    lc = LightCurve(np.arange(n) * dt, sig)
    out = make_whitener(lc)(lc)
    assert out.signal.mean() == pytest.approx(sig.mean(), rel=1e-6)
    f = np.fft.rfftfreq(n, dt)
    p_in, p_out = np.abs(np.fft.rfft(sig - sig.mean())) ** 2, np.abs(np.fft.rfft(out.signal - out.signal.mean())) ** 2
    low, high = (f > 0.05) & (f < 0.5), f > 100
    assert p_out[low].mean() / p_out[high].mean() < 0.1 * p_in[low].mean() / p_in[high].mean()


def test_threshold_for_far():
    snr = np.array([3.2, 3.5, 4.0, 6.0, 8.0])
    assert threshold_for_far(snr, 1.0, 5.0, 3.0) == 3.0           # all 5 allowed
    assert threshold_for_far(snr, 1.0, 0.0, 3.0) > 8.0            # none allowed
    thr = threshold_for_far(snr, 2.0, 1.0, 3.0)                   # 2 allowed in 2 h
    assert (snr >= thr).sum() == 2


def test_veto_mask():
    m = Method("v", lambda lc: lc, veto_times=np.array([10.0]), veto_pad_s=1.0)
    assert list(_vetoed(np.array([8.0, 9.5, 10.5, 12.0]), m)) == [False, True, True, False]


# ───────────────────────────── array injections ─────────────────────────────

from kbo_occultation.analysis import min_size as _ms
from kbo_occultation.analysis.noise_removal import (ComparisonConfig, injection_efficiency_array,
                                                    shadow_offsets_s)
from kbo_occultation.analysis.mimic import MimicConfig, bin_lightcurve

_DT = 5e-4


def _array_setup(depth=0.05, n=120_000, seed=0, sigma=0.01):
    rng = np.random.default_rng(seed)
    tt = np.arange(-0.4, 0.4, 1e-4)
    ti = 1.0 - depth * np.exp(-(tt / 0.05) ** 2)
    v = 25000.0
    time = 1.7e9 + np.arange(n) * _DT
    lcs = {ch: LightCurve(time, 1.0 + sigma * rng.standard_normal(n)) for ch in "ABC"}
    methods = {ch: Method("raw", lambda lc: lc) for ch in lcs}
    patterns = {ch: {100.0: (tt * v, ti)} for ch in lcs}
    templates_m = {ch: {100.0: (tt, ti)} for ch in lcs}
    cfg = ComparisonConfig(mimic=MimicConfig(bin_samples=10), n_trials=6, segment_s=20.0, seed=1)
    return lcs, methods, patterns, templates_m, tt, ti, v, cfg


def test_array_injection_peak_snr_matches_expected_snr_on_white_noise():
    lcs, methods, patterns, templates_m, tt, ti, v, cfg = _array_setup()
    sf = {"A": 0.9, "B": 0.6, "C": 0.8}
    trials, cands = injection_efficiency_array(
        lcs, methods, patterns, templates_m, 0.8, v, sf, [8.0, 52.0], cfg, np.random.default_rng(2),
        [100.0], direction_enu=np.array([1.0, 0.0, 0.0]))
    assert len(trials) == cfg.n_trials * 3
    t_bin = bin_lightcurve(LightCurve(1.7e9 + np.arange(len(tt)) * 1e-4 + 0.0, np.ones(len(tt))), 1)  # noqa: F841
    # expected: binned template residual (star-diluted) over the binned white noise
    dt_b = 10 * _DT
    grid = np.arange(tt.min(), tt.max(), dt_b)
    res = np.interp(grid, tt, ti) - 1.0
    freq = np.linspace(0.0, 0.5 / dt_b, 4097)
    for ch in "ABC":
        sigma_b = 0.01 / np.sqrt(10)
        psd = np.full_like(freq, 2.0 * sigma_b ** 2 * dt_b)
        expect = _ms.expected_snr(sf[ch] * res, freq, psd, 1.0 / dt_b)
        got = np.median([r["peak_snr"] for r in trials if r["channel"] == ch])
        assert got == pytest.approx(expect, rel=0.15)


def test_array_injection_good_window_of_60s_is_finite_and_zero_depth_is_not_detected():
    lcs, methods, patterns, templates_m, tt, ti, v, cfg = _array_setup()
    patterns = {ch: {100.0: (tt * v, np.ones_like(ti))} for ch in lcs}          # nothing injected, dip bank kept
    trials, cands = injection_efficiency_array(
        lcs, methods, patterns, templates_m, 0.8, v, {c: 1.0 for c in lcs}, [8.0, 52.0], cfg,
        np.random.default_rng(3), [100.0])
    assert len(trials) == cfg.n_trials * 3
    assert np.all(np.isfinite([r["peak_snr"] for r in trials]))
    assert max(r["peak_snr"] for r in trials) < 6.0


def test_array_injection_too_short_good_window_returns_nothing():
    lcs, methods, patterns, templates_m, tt, ti, v, cfg = _array_setup()
    trials, cands = injection_efficiency_array(
        lcs, methods, patterns, templates_m, 0.8, v, {c: 1.0 for c in lcs}, [20.0, 22.0], cfg,
        np.random.default_rng(3), [100.0])
    assert trials == [] and cands == []


def test_shadow_offsets_follow_baseline_projection():
    from kbo_occultation.config import ORM_ARRAY
    off = shadow_offsets_s(np.array([1.0, 0.0, 0.0]), ["A", "B", "C"], 25000.0)
    assert off["A"] == 0.0
    assert off["B"] == pytest.approx(ORM_ARRAY.baseline_enu_m("A", "B")[0] / 25000.0)
    assert off["C"] == pytest.approx(ORM_ARRAY.baseline_enu_m("A", "C")[0] / 25000.0)
    iso = shadow_offsets_s(None, ["A", "B", "C"], 25000.0, np.random.default_rng(0))
    assert all(abs(x) <= 200.0 / 25000.0 for x in iso.values())


def test_search_lightcurve_peak_near_sits_at_the_injected_centre():
    from kbo_occultation.analysis.mimic import search_lightcurve
    lcs, methods, patterns, templates_m, tt, ti, v, cfg = _array_setup(depth=0.2, n=60_000)
    lc = lcs["A"]
    t0 = float(lc.time[0] + 15.0)
    inj = inject_occultation(lc.time, lc.signal, tt * v, ti, v, t0_s=t0)
    binned = bin_lightcurve(LightCurve(lc.time, inj.injected_signal), 10)
    rows, _, calib = search_lightcurve(binned, templates_m["A"], 3.0, 0.8, None, peak_near=(t0, 0.05))
    assert abs(calib["peak_time"] - t0) <= 0.01
    assert calib["peak_snr"] == pytest.approx(max(r["snr"] for r in rows), rel=1e-9)
    assert calib["peak_radius_m"] == 100.0
    # a window with no lag gives nan, and no window leaves calib untouched
    assert np.isnan(search_lightcurve(binned, templates_m["A"], 3.0, 0.8, None, peak_near=(t0 + 1e5, 0.05))[2]["peak_snr"])
    assert "peak_snr" not in search_lightcurve(binned, templates_m["A"], 3.0, 0.8, None)[2]


def test_array_injection_deep_event_is_always_detected_in_coincidence_and_shifted_by_the_offsets():
    from kbo_occultation.config import ORM_ARRAY
    lcs, methods, patterns, templates_m, tt, ti, v, cfg = _array_setup(depth=0.3)
    v_slow = 2000.0                                       # slow shadow: offsets of tens of ms, above the 5 ms bin
    patterns = {ch: {100.0: (tt * v_slow, ti)} for ch in lcs}
    trials, cands = injection_efficiency_array(
        lcs, methods, patterns, templates_m, 0.8, v_slow, {c: 1.0 for c in lcs}, [8.0, 52.0], cfg,
        np.random.default_rng(4), [100.0], direction_enu=np.array([1.0, 0.0, 0.0]))
    t = pd.DataFrame(trials)
    c = pd.DataFrame(cands)
    for df in (t, c):
        df["night"], df["observation"], df["correction"], df["sf_mode"] = "n", "o", "raw", "measured"
    # the pattern really sits at t0 + offset: the peak is at the shifted time, not at the common t0
    assert (t["peak_dt_s"].abs() <= 0.01).all()
    for k, g in t.groupby("trial"):
        g = g.set_index("channel")
        assert g.loc["B", "t0_s"] - g.loc["A", "t0_s"] == pytest.approx(ORM_ARRAY.baseline_enu_m("A", "B")[0] / v_slow, abs=1e-6)
        assert g.loc["C", "t0_s"] - g.loc["A", "t0_s"] == pytest.approx(ORM_ARRAY.baseline_enu_m("A", "C")[0] / v_slow, abs=1e-6)
    thr = pd.DataFrame([{"far_per_hr": 1.0, "config": "2of3", "k": 2, "threshold": 6.0},
                        {"far_per_hr": 1.0, "config": "3of3", "k": 3, "threshold": 6.0}])
    eff = _ms.measured_efficiency(t, c, thr, tolerance_s=0.4)
    assert (eff["eff"] == 1.0).all() and (eff["n"] == cfg.n_trials).all()
