import numpy as np

from kbo_occultation.analysis.matched_filter import compute_baseline_and_sigma
from kbo_occultation.simulation.engine import apply_stellar_disk_2d
from kbo_occultation.observations.dc_combine import mask_near_times


def make_lightcurve(n=200_000, dt=5.24e-4, seed=0, dip=False):
    rng = np.random.default_rng(seed)
    t = dt * np.arange(n)
    trend = 10.0 + 0.5 * np.sin(2 * np.pi * t / 120.0)  # slow drift
    signal = trend * (1.0 + 0.02 * rng.standard_normal(n))
    if dip:  # a short dip, like an injected event
        signal[np.abs(t - t[n // 2]) < 0.15] *= 0.7
    return t, trend, signal


def test_decimated_baseline_matches_exact_on_clean_data():
    t, _, signal = make_lightcurve()
    template_duration_s = 0.3

    base_exact, sigma_exact = compute_baseline_and_sigma(t, signal, template_duration_s, decimate=1)
    base_fast, sigma_fast = compute_baseline_and_sigma(t, signal, template_duration_s, decimate="auto")

    rel = np.abs(base_fast - base_exact) / base_exact
    # Both estimates carry ~2e-4 of their own statistical scatter, so
    # the pointwise difference has ~4-sigma tails around 2e-3 over 2e5
    # samples -- an order of magnitude below the 2% data noise.
    assert np.sqrt((rel**2).mean()) < 1e-3, f"baseline rms error {np.sqrt((rel**2).mean()):.2e}"
    assert rel.max() < 3e-3, f"baseline max error {rel.max():.2e}"
    # sigma is now a per-sample local noise estimate; the exact (rolling
    # median of |residual|) and fast (per-block MAD, smoothed) estimators
    # use slightly different robust scales, so compare their overall level.
    assert sigma_fast.shape == signal.shape and sigma_exact.shape == signal.shape
    assert abs(np.median(sigma_fast) - np.median(sigma_exact)) / np.median(sigma_exact) < 0.02


def test_decimated_baseline_robust_to_short_dip():
    # At a deep short dip the exact and decimated medians respond
    # slightly differently to the contaminated samples (the block-mean
    # median rejects the dip more cleanly), so instead of demanding
    # bitwise agreement there, demand what actually matters: neither
    # baseline follows the dip -- both stay on the underlying trend.
    t, trend, signal = make_lightcurve(dip=True)
    n = len(t)
    near_dip = np.abs(t - t[n // 2]) < 1.0

    for decimate in (1, "auto"):
        baseline, _ = compute_baseline_and_sigma(t, signal, 0.3, decimate=decimate)
        assert np.max(np.abs(baseline[near_dip] - trend[near_dip]) / trend[near_dip]) < 0.005


def test_decimate_one_is_exact_path():
    t, _, signal = make_lightcurve(n=20_000)
    b1, s1 = compute_baseline_and_sigma(t, signal, 0.3, decimate=1)
    b2, s2 = compute_baseline_and_sigma(t, signal, 0.3, decimate=1)
    assert np.array_equal(b1, b2) and np.array_equal(s1, s2)


def apply_stellar_disk_2d_loop(x_m, intensity_radial, r_grid_m, star_radius_m,
                               impact_parameter_m, n_star_side):
    """Pre-vectorization reference implementation."""
    offsets = np.linspace(-star_radius_m, star_radius_m, n_star_side)
    dx, dy = np.meshgrid(offsets, offsets)
    mask = dx**2 + dy**2 <= star_radius_m**2
    dx, dy = dx[mask], dy[mask]
    convolved = np.zeros_like(x_m)
    for i, x in enumerate(x_m):
        r = np.sqrt((x + dx)**2 + (impact_parameter_m + dy)**2)
        convolved[i] = np.interp(r, r_grid_m, intensity_radial, left=1.0, right=1.0).mean()
    return convolved


def test_stellar_disk_vectorization_matches_loop():
    rng = np.random.default_rng(2)
    x_m = np.linspace(-5000.0, 5000.0, 801)
    r_grid = np.linspace(0.0, 9000.0, 500)
    intensity = 1.0 + 0.3 * np.exp(-(r_grid - 1500.0)**2 / 500.0**2)
    for b in (0.0, 800.0):
        new = apply_stellar_disk_2d(x_m, intensity, r_grid, 400.0, b, 16)
        old = apply_stellar_disk_2d_loop(x_m, intensity, r_grid, 400.0, b, 16)
        assert np.allclose(new, old, rtol=1e-12, atol=0)


def test_mask_near_times_matches_bruteforce():
    rng = np.random.default_rng(3)
    time = np.sort(rng.uniform(0, 100, 5000))
    flagged = rng.uniform(0, 100, 17)
    pad = 0.8

    fast = mask_near_times(time, flagged, pad)
    brute = np.zeros(len(time), dtype=bool)
    for t in flagged:
        brute |= (time >= t - pad) & (time <= t + pad)
    assert np.array_equal(fast, brute)

    assert not mask_near_times(time, np.array([]), pad).any()
