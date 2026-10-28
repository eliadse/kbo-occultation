import numpy as np
import pytest

from kbo_occultation.analysis.snr import snr_report
from kbo_occultation.observations.log import load_observation_log
from kbo_occultation.observations.prepare import build_caches

from test_photometry_loader import write_bin

STEM = "Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_{}_10002_20002_{}"
N = 60_000                       # ~31 s at the standard cadence
STAR_VAR, NSB_VAR, REL_SIGMA = 10.0, 1.0, 0.01

LOG = """
site: {{extinction_k: {{B: 0.23, V: 0.13}}}}
nights:
  "2026-01-01":
    date: "2026-01-01"
    channels: [A]
    dark_run: dark
    observations:
      HD-1:
        file: {star}.npz
        type: star
        ra_deg: 90.0
        dec_deg: 30.0
        mag_B: 7.0
        mag_V: 7.0
        good_time: null
      dark:
        file: {dark}.npz
        type: dark
        good_time: null
"""


@pytest.fixture
def white_noise_log(tmp_path, monkeypatch):
    monkeypatch.delenv("KBO_OBSERVATIONS_DIR", raising=False)
    rng = np.random.default_rng(3)
    star, dark = STEM.format("HD-1", "20260101T220000"), STEM.format("dark", "20260101T230000")
    for stem, level in ((star, STAR_VAR), (dark, NSB_VAR)):
        var = rng.normal(level, REL_SIGMA * level, N)
        write_bin(tmp_path / f"{stem}.bin", np.sqrt(var).astype(np.float32))
    (tmp_path / "observation_log.yaml").write_text(LOG.format(star=star, dark=dark))
    log = load_observation_log(str(tmp_path / "observation_log.yaml"))
    build_caches(log.nights["2026-01-01"], log=lambda *_: None)
    return log


def test_snr_report_on_white_noise(white_noise_log):
    res = snr_report(white_noise_log, n_bins=(10, 50), detrend_window_s=2.0)
    assert list(res.index) == [("2026-01-01", "HD-1", "MAGIC-2", 10), ("2026-01-01", "HD-1", "MAGIC-2", 50)]
    row = res.loc[("2026-01-01", "HD-1", "MAGIC-2", 10)]
    assert row["nsb_source"] == "dark"

    # white noise: the bin-to-bin SNR is (S - NSB) / (sigma / sqrt(n_bin))
    expected = (STAR_VAR - NSB_VAR) / (REL_SIGMA * STAR_VAR / np.sqrt(10))
    assert row["detr"] == pytest.approx(expected, rel=0.1)
    assert row["mean_fast"] == pytest.approx(expected, rel=0.15)

    # no DC report -> no DC-based columns; theory present at the computed airmass
    assert "mean_comb" not in res or np.isnan(row["mean_comb"])
    assert row["airmass"] >= 1.0 and row["mirror_throughput"] == 1.0
    assert row["theo_obs_B"] > 0 and row["ratio_detr_B"] == pytest.approx(row["detr"] / row["theo_obs_B"])


def test_snr_report_has_inband_noise_columns(white_noise_log):
    res = snr_report(white_noise_log, n_bins=(10,), detrend_window_s=2.0)
    row = res.iloc[0]
    for col in ("inband_rms", "mf_sigma", "welch_short"):
        assert col in res.columns
    # white noise of per-sample rel. sigma REL_SIGMA: a filter sees the full sigma
    # at the variance-flux level, the band holds only a fraction of it
    assert row["mf_sigma"] == pytest.approx(REL_SIGMA, rel=0.3)
    assert 0 < row["inband_rms"] < REL_SIGMA
