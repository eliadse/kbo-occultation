import os

import numpy as np
import pytest

from kbo_occultation.observations.cc_reports import (read_cc_reports, summarize_conditions,
                                         save_conditions, load_conditions, pointing_conditions)
from kbo_occultation.observations.log import load_observation_log

FIXTURE = os.path.join(os.path.dirname(__file__), "data", "cc_excerpt.rep")


def test_read_cc_reports_fields():
    rep = read_cc_reports(FIXTURE)
    cc = rep["cc"]
    # 2026-03-08 22:15 line: HD-40626 at Zd 31.73, Az -39.50; the 2025 "Park"
    # line has 00.00 pointing, which must come out as missing, not zenith.
    hd = cc[cc["source"] == "HD-40626"].iloc[0]
    assert hd["zd_deg"] == pytest.approx(31.73)
    assert hd["az_deg"] == pytest.approx(-39.50)
    assert hd["target_ra_deg"] == pytest.approx(90.70)
    assert hd["humidity_pct"] == pytest.approx(24.12)
    assert hd["moon_flag"] == "No_Moon"
    assert np.isnan(cc[cc["source"] == "Park"].iloc[0]["zd_deg"])

    w = rep["weather"].iloc[0]
    assert w["temperature_C"] == pytest.approx(4.57)
    assert w["humidity_pct"] == pytest.approx(23.59)
    assert w["wind_avg_kmh"] == pytest.approx(19.72)

    p = rep["pyro"].iloc[0]
    assert p["sky_temperature_K"] == pytest.approx(208.11)
    assert p["cloudiness_pct"] == pytest.approx(12.54)

    a = rep["amc"].iloc[0]
    assert [a[f"amc_panels_{k}"] for k in range(4)] == [233, 7, 0, 5]

    s = rep["starg"].iloc[0]
    assert s["starg_n_stars"] == 20
    assert s["starg_sky_brightness"] == pytest.approx(88.9)

    assert rep["lidar"].iloc[0]["lidar_error_code"] == 4


def test_summary_and_cache_roundtrip(tmp_path):
    rep = read_cc_reports(FIXTURE)
    t = rep["cc"]["time"].max()
    save_conditions(rep, str(tmp_path), "m2")
    loaded = load_conditions(str(tmp_path))["m2"]
    summary = summarize_conditions(loaded, t - 60, t + 60)
    assert summary["cc_source"] == "HD-40626"
    assert summary["cc_zd_mean_deg"] == pytest.approx(31.73)
    assert summary["cloudiness_pct"] == pytest.approx(12.54)
    assert summary["lidar_operating"] is False


def test_pointing_matches_cc_report():
    # HD-40626 at 2026-03-08T22:15:07 UTC: the CC report gives Zd 31.73 deg.
    from astropy.time import Time
    obs = load_observation_log().nights["2026-03-08"].observations["HD-40626"]
    t = Time("2026-03-08T22:15:07", scale="utc").unix
    p = pointing_conditions(obs.ra_deg, obs.dec_deg, t, t + 1, n_points=2)
    assert p["zd_mean_deg"] == pytest.approx(31.73, abs=0.2)
    assert p["moon_alt_deg"] < 0  # CC flag: No_Moon
