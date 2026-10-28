import numpy as np
import pytest

from kbo_occultation.analysis.matched_filter import (
    calibrate_chi2_threshold, empirical_filter_sigma, sliding_matched_filter_snr)

DT = 0.01
N = 60000  # 600 s


def _template():
    t = np.arange(-0.5, 0.5, DT)
    return t, 1.0 - 0.3 * np.exp(-0.5 * (t / 0.1) ** 2) * np.cos(6 * t) ** 2


def _red_noise(rng, amp=0.02, f_lo=0.05):
    """1/f-ish noise between f_lo and Nyquist, unit-scaled to rms `amp`."""
    f = np.fft.rfftfreq(N, DT)
    spec = np.where(f >= f_lo, 1.0 / np.sqrt(np.maximum(f, f_lo)), 0.0)
    ph = np.exp(2j * np.pi * rng.random(len(f)))
    x = np.fft.irfft(spec * ph, n=N)
    return amp * x / x.std()


def _snr_std(noise, normalization, **kw):
    t = np.arange(N) * DT
    tt, ti = _template()
    res = sliding_matched_filter_snr(t, 1.0 + noise, tt, ti, normalization=normalization, **kw)
    return np.std(res.snr[200:-200]), res


def test_white_noise_all_normalizations_agree():
    rng = np.random.default_rng(1)
    noise = 0.01 * rng.standard_normal(N)
    stds = {n: _snr_std(noise, n)[0] for n in ("local_block", "global", "empirical")}
    for v in stds.values():
        assert v == pytest.approx(1.0, rel=0.12), stds
    assert max(stds.values()) / min(stds.values()) < 1.12


def test_red_noise_empirical_snr_has_unit_sigma():
    rng = np.random.default_rng(2)
    noise = _red_noise(rng)
    emp, res = _snr_std(noise, "empirical")
    loc, _ = _snr_std(noise, "local_block")
    assert emp == pytest.approx(1.0, rel=0.25)
    # the white-noise normalisation is optimistic on red noise
    assert loc > 1.5 * emp
    assert res.amp_sigma is not None and res.normalization == "empirical"


def test_empirical_sigma_ignores_an_event():
    rng = np.random.default_rng(3)
    x = rng.standard_normal(20000)
    clean = empirical_filter_sigma(x, 2000)
    x[10000:10100] += 30.0
    dirty = empirical_filter_sigma(x, 2000)
    assert np.allclose(clean, dirty, rtol=0.05)


def test_unknown_normalization_raises():
    t = np.arange(N) * DT
    tt, ti = _template()
    with pytest.raises(ValueError):
        sliding_matched_filter_snr(t, np.ones(N), tt, ti, normalization="bogus")


def test_chi2_calibration_near_one_on_white_noise():
    rng = np.random.default_rng(4)
    t = np.arange(N) * DT
    tt, ti = _template()
    res = sliding_matched_filter_snr(t, 1.0 + 0.01 * rng.standard_normal(N), tt, ti)
    med, q99 = calibrate_chi2_threshold(res, weighting="empirical", n_samples=100)
    assert med == pytest.approx(1.0, rel=0.15)
    assert q99 > med


def test_spectral_matched_filter_sigma_matches_empirical():
    """Part 3.1 target: PSD-based sigma within 10 % of the measured filter-output scatter."""
    from kbo_occultation.analysis.spectra import matched_filter_sigma, segmented_welch

    rng = np.random.default_rng(5)
    noise = _red_noise(rng)
    t = np.arange(N) * DT
    tt, ti = _template()
    res = sliding_matched_filter_snr(t, 1.0 + noise, tt, ti, normalization="empirical")
    f, p = segmented_welch(t, res.residual_data, 1 / DT, 2 ** 13)
    sigma_psd = matched_filter_sigma(f, p, res.residual_template, 1 / DT)
    amp = res.snr * res.amp_sigma                 # filter output, flux units
    sigma_emp = float(np.std(amp[200:-200]))
    assert sigma_psd == pytest.approx(sigma_emp, rel=0.10)


def test_sigma_path_gap_on_red_noise():
    """decimate="auto" vs 1: local_block sigma depends on the path, empirical does not."""
    from kbo_occultation.analysis.matched_filter import compute_baseline_and_sigma

    rng = np.random.default_rng(6)
    noise = _red_noise(rng)
    t = np.arange(N) * DT
    tt, ti = _template()
    std = {}
    sigma_med = {}
    for dec in ("auto", 1):
        b, s = compute_baseline_and_sigma(t, 1.0 + noise, 1.0, decimate=dec)
        sigma_med[dec] = np.median(s)
        for norm in ("local_block", "empirical"):
            r = sliding_matched_filter_snr(t, 1.0 + noise, tt, ti, baseline=b, sigma=s, normalization=norm)
            std[(dec, norm)] = np.std(r.snr[200:-200])
    # per-block MAD of short blocks sees only the high-frequency part of red noise
    assert sigma_med["auto"] < 0.9 * sigma_med[1]
    assert std[("auto", "local_block")] > 1.3 * std[(1, "local_block")] > 2.0
    # the empirical normalisation does not care which sigma path was used
    assert std[("auto", "empirical")] == pytest.approx(std[(1, "empirical")], rel=0.05)
    assert std[("auto", "empirical")] == pytest.approx(1.0, rel=0.2)
