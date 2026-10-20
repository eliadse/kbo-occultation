import numpy as np

from kbo_occultation.analysis.injection import inject_occultation


def make_template():
    x_m = np.linspace(-3000.0, 3000.0, 500)
    intensity = 1.0 - 0.4 * np.exp(-x_m**2 / 500.0**2)
    return x_m, intensity


def test_float32_signal_stays_float32():
    x_m, intensity = make_template()
    t = 5.24e-4 * np.arange(50_000)
    signal32 = np.full(len(t), 10.0, dtype=np.float32)

    inj = inject_occultation(t, signal32, x_m, intensity, 25_000.0, t0_s=t[len(t) // 2])
    assert inj.injected_signal.dtype == np.float32
    assert inj.template.dtype == np.float32


def test_float32_matches_float64_path():
    x_m, intensity = make_template()
    t = 5.24e-4 * np.arange(50_000)
    rng = np.random.default_rng(0)
    base = 10.0 + rng.standard_normal(len(t))
    t0 = t[len(t) // 2]

    inj64 = inject_occultation(t, base.astype(np.float64), x_m, intensity, 25_000.0, t0_s=t0)
    inj32 = inject_occultation(t, base.astype(np.float32), x_m, intensity, 25_000.0, t0_s=t0)
    assert np.allclose(inj32.injected_signal, inj64.injected_signal, rtol=1e-5)
    assert inj32.window_s == inj64.window_s


def test_injection_is_multiplicative_and_windowed():
    x_m, intensity = make_template()
    t = 5.24e-4 * np.arange(50_000)
    signal = np.full(len(t), 4.0)
    t0 = t[len(t) // 2]

    inj = inject_occultation(t, signal, x_m, intensity, 25_000.0, t0_s=t0)
    # well outside the window the template is 1 (the window edge itself
    # sits at the |template - 1| > 1e-3 tolerance, so leave some margin
    # for the template's own tails)
    lo, hi = inj.window_s
    outside = (t < lo - 0.05) | (t > hi + 0.05)
    assert np.allclose(inj.injected_signal[outside], 4.0, rtol=1e-4)
    # at closest approach the dip depth matches the template minimum
    assert np.isclose(inj.injected_signal.min(), 4.0 * intensity.min(), rtol=1e-3)
