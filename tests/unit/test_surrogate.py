import numpy as np
import pytest

from kbo_occultation.analysis.spectra import segmented_welch, inband_rms, matched_filter_sigma
from kbo_occultation.analysis.surrogates import (
    phase_randomize_array, phase_randomized_lightcurve, surrogate_lightcurves)
from kbo_occultation.observations.photometry import LightCurve

FS = 1000.0


def _red(rng, n):
    f = np.fft.rfftfreq(n, 1 / FS)
    spec = 1.0 / np.sqrt(np.maximum(f, 0.2))
    x = np.fft.irfft(spec * np.exp(2j * np.pi * rng.random(len(f))) * np.sqrt(n), n=n)
    return 1.0 + 0.02 * x / x.std()


def _gapped(rng, n=2 ** 18):
    t = np.arange(n) / FS
    sig = _red(rng, n)
    keep = np.ones(n, bool)
    keep[60000:61500] = False        # a hole, like despike(fill="remove")
    keep[150000:151000] = False
    return LightCurve(t[keep], sig[keep])


def test_amplitude_spectrum_and_mean_preserved():
    rng = np.random.default_rng(0)
    x = _red(rng, 2 ** 16)
    s = phase_randomize_array(x, rng)
    assert s.mean() == pytest.approx(x.mean())
    assert np.allclose(np.abs(np.fft.rfft(s - s.mean())), np.abs(np.fft.rfft(x - x.mean())))
    assert not np.allclose(s, x)


def test_psd_preserved_within_welch_error():
    rng = np.random.default_rng(1)
    lc = _gapped(rng)
    sur = phase_randomized_lightcurve(lc, rng)
    f, p0 = segmented_welch(lc.time, lc.signal / lc.signal.mean() - 1, FS, 2 ** 12)
    _, p1 = segmented_welch(sur.time, sur.signal / sur.signal.mean() - 1, FS, 2 ** 12)
    band = (f > 0.5) & (f < 200)
    ratio = p1[band] / p0[band]
    assert abs(np.median(ratio) - 1) < 0.05
    assert abs(np.exp(np.mean(np.log(ratio))) - 1) < 0.1


def test_gap_mask_and_runs_kept():
    rng = np.random.default_rng(2)
    lc = _gapped(rng)
    sur = phase_randomized_lightcurve(lc, rng)
    assert np.array_equal(sur.time, lc.time)
    gaps = np.where(np.diff(lc.time) > 1.5 / FS)[0]
    edges = np.concatenate([[0], gaps + 1, [len(lc.time)]])
    for a, b in zip(edges[:-1], edges[1:]):
        assert sur.signal[a:b].mean() == pytest.approx(lc.signal[a:b].mean())


def test_surrogates_are_independent_and_reproducible():
    lc = _gapped(np.random.default_rng(3), 2 ** 15)
    a = [s.signal for s in surrogate_lightcurves(lc, 2, seed=5)]
    b = [s.signal for s in surrogate_lightcurves(lc, 2, seed=5)]
    assert np.array_equal(a[0], b[0]) and not np.array_equal(a[0], a[1])


def test_matched_filter_sigma_white_and_inband():
    rng = np.random.default_rng(4)
    x = 0.01 * rng.standard_normal(2 ** 18)
    t = np.arange(len(x)) / FS
    f, p = segmented_welch(t, x, FS, 2 ** 14)
    r = -0.3 * np.exp(-0.5 * (np.arange(-300, 300) / 50.0) ** 2)
    assert matched_filter_sigma(f, p, r, FS) == pytest.approx(0.01, rel=0.05)
    lc = LightCurve(t, 1.0 + x)
    # white noise of sigma s has rms s*sqrt((f2-f1)/(fs/2)) in a band
    assert inband_rms(lc, (10.0, 60.0), 2 ** 14) == pytest.approx(0.01 * np.sqrt(50 / 500), rel=0.08)
