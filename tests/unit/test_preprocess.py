import numpy as np

from kbo_occultation.observations.photometry import (LightCurve, preprocess_observation)

from test_photometry_loader import write_bin, DT_S


def test_preprocess_roundtrip_and_cleaning(tmp_path):
    # A flat channel with a single-sample positive spike (a cosmic-hit style
    # glitch) that the outlier correction should remove.
    rng = np.random.default_rng(0)
    std = rng.normal(3.0, 0.05, 4000).astype(np.float32)
    spike_i = 2000
    std[spike_i] = 100.0  # huge positive excursion -> variance spike
    fn = write_bin(tmp_path / "obs.bin", std)

    out = preprocess_observation(str(fn), str(tmp_path / "obs.npz"))
    assert out == str(tmp_path / "obs.npz")

    # raw variance is exactly std**2 and preserves the spike
    raw = LightCurve.from_preprocessed(out, "A", cleaned=False)
    assert np.allclose(raw.signal, std.astype(float) ** 2, rtol=1e-5)
    assert raw.signal[spike_i] == max(raw.signal)

    # reconstructed grid at the standard cadence
    assert np.allclose(np.diff(raw.time), DT_S)

    # cleaned variance removes the isolated spike (default positive polarity)
    clean = LightCurve.from_preprocessed(out, "A", cleaned=True)
    assert clean.signal[spike_i] < raw.signal[spike_i]
    assert clean.meta["cleaned"] is True


def test_from_preprocessed_all_matches_binary(tmp_path):
    rng = np.random.default_rng(1)
    std_a = rng.uniform(1.0, 2.0, 500).astype(np.float32)
    std_b = rng.uniform(1.0, 2.0, 500).astype(np.float32)
    std_c = rng.uniform(1.0, 2.0, 500).astype(np.float32)
    fn = write_bin(tmp_path / "obs.bin", std_a, std_b, std_c)

    out = preprocess_observation(str(fn), str(tmp_path / "obs.npz"))

    from_bin = LightCurve.from_stat_binary_all(str(fn))
    from_npz = LightCurve.from_preprocessed_all(out)  # cleaned=False by default
    for ch in ["A", "B", "C"]:
        # float32 storage in the cache -> compare with a tolerance
        assert np.allclose(from_npz[ch].signal, from_bin[ch].signal, rtol=1e-5)
        assert np.allclose(from_npz[ch].time, from_bin[ch].time)
