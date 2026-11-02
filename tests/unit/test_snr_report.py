import numpy as np
import pytest

from kbo_occultation.analysis.snr import load_lightcurves, nsb_table, snr_report
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


def test_nsb_table_units_and_star_peaks(white_noise_log):
    fast = nsb_table(white_noise_log, units="fast").iloc[0]
    assert fast["units"] == "fast" and fast["nsb_source"] == "dark"
    assert fast["star_mean"] == pytest.approx(STAR_VAR, rel=0.01)
    # on white noise the "last" peak stops at the first noisy local maximum from the top (~3 sigma)
    assert fast["star_peak_fast"] == pytest.approx(STAR_VAR, rel=0.05)
    assert fast["nsb_default"] == fast["nsb_default_fast"] == pytest.approx(NSB_VAR, rel=0.01)

    # no DC report: no dc_scale, so DC-unit levels are NaN but the fast default survives
    dc = nsb_table(white_noise_log).iloc[0]
    assert dc["units"] == "dc"
    assert np.isnan(dc["star_mean"]) and np.isnan(dc["nsb_default"]) and np.isnan(dc["star_peak_comb"])
    assert dc["nsb_default_fast"] == pytest.approx(NSB_VAR, rel=0.01)

    with pytest.raises(ValueError):
        nsb_table(white_noise_log, units="volts")


DC_SCALE = 4.0


@pytest.fixture
def synthetic_dc(white_noise_log, monkeypatch):
    """A 1 Hz DC report: on-target pixel follows the runs' levels / DC_SCALE, background pixels sit at the NSB."""
    from kbo_occultation.analysis import snr, snr_measurement
    from kbo_occultation.observations import dc_combine

    from kbo_occultation.observations.photometry import LightCurve

    lcs = load_lightcurves(white_noise_log)
    dark = lcs[("2026-01-01", "dark", "A")]          # the test runs share time stamps: move the dark run 1 h later
    lcs[("2026-01-01", "dark", "A")] = LightCurve(dark.time + 3600.0, dark.signal, meta=dict(dark.meta))
    star_lc = lcs[("2026-01-01", "HD-1", "A")]
    t0, t1 = min(lc.time[0] for lc in lcs.values()) - 60, max(lc.time[-1] for lc in lcs.values()) + 60
    t = np.arange(t0, t1, 1.0)
    rng = np.random.default_rng(5)

    def fake_load_dc_series(path, telescope="m2", pixel=0):
        on_star = (t >= star_lc.time[0]) & (t <= star_lc.time[-1])
        level = np.where(on_star & (pixel == 0), STAR_VAR, NSB_VAR) / DC_SCALE
        return t, level * (1 + REL_SIGMA * rng.standard_normal(t.size))

    monkeypatch.setattr(dc_combine, "load_dc_series", fake_load_dc_series)
    monkeypatch.setattr(snr_measurement, "load_dc_series", fake_load_dc_series)
    monkeypatch.setattr(snr, "_has_dc", lambda night: True)
    return lcs


def test_nsb_table_dark_comb_peak(white_noise_log, synthetic_dc):
    row = nsb_table(white_noise_log, lcs=synthetic_dc).iloc[0]
    assert row["nsb_source"] == "dark_comb_peak" and row["units"] == "dc"
    assert row["dc_scale"] == pytest.approx(DC_SCALE, rel=0.02)
    assert row["nsb_default"] == row["nsb_dark_comb_peak"] == pytest.approx(NSB_VAR / DC_SCALE, rel=0.01)
    assert row["star_peak_comb"] == pytest.approx(STAR_VAR / DC_SCALE, rel=0.05)
    assert row["star_peak_fast"] == pytest.approx(STAR_VAR / DC_SCALE, rel=0.05)       # converted to DC units
    assert row["nsb_default_fast"] == pytest.approx(NSB_VAR, rel=0.05)

    fast = nsb_table(white_noise_log, lcs=synthetic_dc, units="fast").iloc[0]
    assert fast["nsb_dark_comb_peak"] == pytest.approx(NSB_VAR, rel=0.02)
    assert fast["nsb_dark"] == pytest.approx(NSB_VAR, rel=0.01)


def test_gaussian_histogram_peak_ignores_dip_tail():
    from kbo_occultation.analysis.snr_measurement import gaussian_histogram_peak

    rng = np.random.default_rng(7)
    core = rng.normal(10.0, 0.2, 9000)
    dips = rng.uniform(6.0, 9.5, 1000)                # 10 % of the samples in a low tail
    mu, sigma, amp = gaussian_histogram_peak(np.concatenate([core, dips]))
    assert mu == pytest.approx(10.0, rel=0.01) and sigma == pytest.approx(0.2, rel=0.15) and amp > 0
    with pytest.warns(UserWarning):
        assert np.isnan(gaussian_histogram_peak([1.0, 2.0])[0])
