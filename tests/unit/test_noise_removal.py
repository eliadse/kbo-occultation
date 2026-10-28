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
    assert new["sigma_floor"] == pytest.approx(np.sqrt(2 * 2 / 2 ** 18), rel=1e-6)   # 0.39 %
    assert new["sigma_total"] > base["sigma_total"]
    assert base["sigma_floor"] == 0.0
    # the floor averages down with the number of frames in a bin, red noise does not
    coarse = tel.noise_budget(6.0, 5 * u.ms, red_noise_rms=2e-3)
    assert coarse["sigma_floor"] < new["sigma_floor"] / 3
    assert coarse["sigma_red"] == pytest.approx(2e-3)
    assert Instrument("LST").noise_budget(6.0, 524.288 * u.us)["sigma_floor"] == 0.0
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
