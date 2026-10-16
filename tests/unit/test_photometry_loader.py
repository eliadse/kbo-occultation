import numpy as np

from kbo_occultation import config
from kbo_occultation.observations.photometry import LightCurve, average_chunks

DATA_FORMAT = [('magic_word', '<Q'),
               ('mean_chA', '<f'), ('mean_chB', '<f'), ('mean_chC', '<f'), ('mean_chD', '<f'),
               ('std_chA', '<f'), ('std_chB', '<f'), ('std_chC', '<f'), ('std_chD', '<f'),
               ('min_chA', '<f'), ('min_chB', '<f'), ('min_chC', '<f'), ('min_chD', '<f'),
               ('max_chA', '<f'), ('max_chB', '<f'), ('max_chC', '<f'), ('max_chD', '<f'),
               ('time_stamp', '<Q')]

DT_S = config.standard_sampling * config.standard_sample_duration / 1e9  # ~0.524 ms


def write_bin(path, std_a, std_b=None, std_c=None):
    n = len(std_a)
    rec = np.zeros(n, dtype=DATA_FORMAT)
    rec['std_chA'] = std_a
    rec['std_chB'] = std_b if std_b is not None else std_a
    rec['std_chC'] = std_c if std_c is not None else std_a
    rec['time_stamp'] = (1e6 * (1000.0 + DT_S * np.arange(n))).astype(np.uint64)
    rec.tofile(str(path))
    return path


def test_from_stat_binary_loads_variance(tmp_path):
    std = np.full(1000, 3.0, dtype=np.float32)
    fn = write_bin(tmp_path / "obs.bin", std)

    # The flux proxy is always the variance (std**2); it is no longer selectable.
    lc = LightCurve.from_stat_binary(str(fn), "A")
    assert np.allclose(lc.signal, 9.0)
    assert lc.meta["channel"] == "A"
    # reconstructed uniform grid at the standard cadence
    assert np.allclose(np.diff(lc.time), DT_S)


def test_loaders_agree(tmp_path):
    rng = np.random.default_rng(1)
    std = rng.uniform(1.0, 2.0, 500).astype(np.float32)
    fn = write_bin(tmp_path / "obs.bin", std)

    single = LightCurve.from_stat_binary(str(fn), "A")
    all_ch = LightCurve.from_stat_binary_all(str(fn))
    assert np.array_equal(single.signal, all_ch["A"].signal)


def test_highpass_uses_effective_fs_after_averaging(tmp_path):
    # A slow sinusoid at 1 Hz plus a constant. With average=4 the cadence
    # is 4x coarser; a 5 Hz cutoff must still remove the 1 Hz component,
    # which only happens if the *effective* fs is used.
    n = 8192
    t = DT_S * np.arange(n)
    std = np.sqrt(10.0 + 5.0 * np.sin(2 * np.pi * 1.0 * t)).astype(np.float32)
    fn = write_bin(tmp_path / "obs.bin", std)

    lc = LightCurve.from_stat_binary(str(fn), "A", average=4, low_freq_cut=5.0)
    # After a correct high-pass, the 1 Hz swing (amplitude 5) is gone:
    # residual fluctuation should be tiny compared to the original.
    assert np.std(lc.signal) < 0.5
    assert len(lc.signal) == n // 4


def test_average_chunks():
    x = np.arange(10.0)
    assert np.array_equal(average_chunks(x, 1), x)
    out = average_chunks(x, 3)
    assert np.allclose(out, [1.0, 4.0, 7.0])  # trailing element dropped
