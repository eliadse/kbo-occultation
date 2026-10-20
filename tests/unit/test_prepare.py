import datetime as dt

import numpy as np
import pandas as pd
import pytest
import yaml

from kbo_occultation.observations.log import load_observation_log
from kbo_occultation.observations.prepare import (detect_channels, find_cc_reports, night_files,
                                                  parse_filename, prepare_night)
from kbo_occultation.observations.simbad import simbad_name

from test_photometry_loader import write_bin

STEM = "Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_{}_10002_20002_{}"

LOG_HEADER = """# Header comment that must survive.
paths:
  cc_reports_root: {cc_root}

site:
  extinction_k: {{B: 0.23, V: 0.13}}

nights:
"""

FAKE_SIMBAD = {
    "HD-1": {"ra_deg": 10.0, "dec_deg": 20.0, "mag_B": 7.1, "mag_V": 7.0, "spectral_type": "B9", "simbad_id": "HD 1"},
    "HD-2": {"ra_deg": 30.0, "dec_deg": 40.0, "mag_B": 8.1, "mag_V": 8.0, "spectral_type": "A0", "simbad_id": "HD 2"},
}


def fake_query_stars(names):
    return {n: FAKE_SIMBAD.get(n) for n in names}


def write_run(directory, name, ts, level_a, seed):
    """Channel A sees light at ``level_a`` (variance); B and C sit at a 0.6 pedestal."""
    rng = np.random.default_rng(seed)
    n = 2000
    std_a = np.sqrt(rng.normal(level_a, 0.01 * level_a, n)).astype(np.float32)
    std_ped = np.sqrt(rng.normal(0.6, 0.002, n)).astype(np.float32)
    write_bin(directory / (STEM.format(name, ts) + ".bin"), std_a, std_ped, std_ped.copy())


@pytest.fixture
def night_dir(tmp_path, monkeypatch):
    monkeypatch.delenv("KBO_OBSERVATIONS_DIR", raising=False)
    monkeypatch.delenv("KBO_CC_REPORTS_ROOT", raising=False)
    write_run(tmp_path, "HD-1", "20260101T220000", 20.0, 1)
    write_run(tmp_path, "dark", "20260101T223000", 5.0, 2)
    write_run(tmp_path, "HD-2", "20260102T010000", 12.0, 3)
    (tmp_path / "cc").mkdir()
    log_path = tmp_path / "observation_log.yaml"
    log_path.write_text(LOG_HEADER.format(cc_root=tmp_path / "cc"))
    return tmp_path, log_path


def prepare(log_path, **kw):
    return prepare_night("2026-01-01", log_path=str(log_path), query_stars=fake_query_stars,
                         log=lambda *_: None, **kw)


def test_parse_filename():
    info = parse_filename(STEM.format("HD-4335", "20251216T212254") + ".bin")
    assert info["name"] == "HD-4335"
    assert info["start_utc"] == dt.datetime(2025, 12, 16, 21, 22, 54)
    assert info["daq"]["sample_duration_ns"] == 2
    assert info["daq"]["range_mV"] == 200 and info["daq"]["impedance_ohm"] == 50
    assert parse_filename("something_else.bin") is None


def test_night_window_spans_midnight(tmp_path):
    for name, ts in [("HD1", "20260819T230000"), ("DarkPatch5", "20260820T005147"),
                     ("HD2", "20260820T130000"), ("HD0", "20260819T110000")]:
        write_bin(tmp_path / (STEM.format(name, ts) + ".bin"), np.ones(10, dtype=np.float32))
    names = [parse_filename(p)["name"] for p in night_files("2026-08-19", str(tmp_path))]
    assert names == ["HD1", "DarkPatch5"]


def test_simbad_name():
    assert simbad_name("HD-4335") == "HD 4335"
    assert simbad_name("HD170263") == "HD 170263"
    assert simbad_name("HIP_12345") == "HIP 12345"


def test_detect_channels():
    # 2025-12-16-like: only A lit, B/C at a constant pedestal
    multi = pd.DataFrame({"A": [61.9, 38.3, 8.1], "B": [0.725, 0.727, 0.725], "C": [0.442, 0.441, 0.442]})
    assert detect_channels(multi) == ["A"]
    # dark patch of 2026-08-19: all lit, C lowest at 2.5
    assert detect_channels(pd.DataFrame({"A": [6.5, 26.6], "B": [13.0, 32.1], "C": [2.5, 6.6]})) == ["A", "B", "C"]
    # single run: level alone decides
    assert detect_channels(pd.DataFrame({"A": [20.0], "B": [0.7], "C": [3.0]})) == ["A", "C"]


def test_find_cc_reports(tmp_path):
    (tmp_path / "2026" / "2026_08_19").mkdir(parents=True)
    for f in ["2026/2026_08_19/CC_M1_2026_08_19_18_55_05.rep", "2026/2026_08_19/CC_M2_2026_08_19_18_55_05.rep",
              "2026/CC_M1_2026_03_08_20_14_52.rep"]:
        (tmp_path / "2026").mkdir(exist_ok=True)
        (tmp_path / f).write_text("")
    assert find_cc_reports("2026-08-19", str(tmp_path)) == {
        "m1": "2026/2026_08_19/CC_M1_2026_08_19_18_55_05.rep",
        "m2": "2026/2026_08_19/CC_M2_2026_08_19_18_55_05.rep"}
    assert find_cc_reports("2026-03-08", str(tmp_path)) == {"m1": "2026/CC_M1_2026_03_08_20_14_52.rep"}
    assert find_cc_reports("2026-01-01", None) == {}


def test_prepare_creates_and_fills_the_night(night_dir):
    _, log_path = night_dir
    problems = prepare(log_path)

    text = log_path.read_text()
    assert text.startswith("# Header comment that must survive.")
    raw = yaml.safe_load(text)["nights"]["2026-01-01"]
    assert raw["channels"] == ["A"]
    assert raw["dark_run"] == "dark"
    assert "m1" not in raw["faulty_mirrors"] and "lst1" not in raw["faulty_mirrors"]   # no light there
    assert raw["observations"]["HD-1"]["ra_deg"] == 10.0
    assert raw["observations"]["HD-2"]["spectral_type"] == "A0"
    assert raw["observations"]["dark"]["type"] == "dark"
    assert "# SIMBAD" in text and "# data" in text

    night = load_observation_log(str(log_path)).nights["2026-01-01"]
    assert all(o.has_cache for o in night.observations.values())
    assert night.observations["HD-2"].mag_B == 8.1
    # what is left is manual (or not derivable without CC reports / DC report)
    assert any("good_time" in p for p in problems)
    assert not any("ra_deg" in p and "HD-" in p for p in problems)


def test_prepare_keeps_manual_values_unless_refresh(night_dir):
    _, log_path = night_dir
    prepare(log_path)
    log_path.write_text(log_path.read_text().replace("ra_deg: 10.0", "ra_deg: 11.0  # corrected by hand"))

    prepare(log_path)
    text = log_path.read_text()
    assert "ra_deg: 11.0  # corrected by hand" in text

    prepare(log_path, dry_run=True, refresh=True)
    assert log_path.read_text() == text            # dry run writes nothing

    prepare(log_path, refresh=True)
    assert yaml.safe_load(log_path.read_text())["nights"]["2026-01-01"]["observations"]["HD-1"]["ra_deg"] == 10.0


def test_prepare_adds_new_runs_to_existing_night(night_dir):
    tmp_path, log_path = night_dir
    prepare(log_path)
    write_run(tmp_path, "HD-3", "20260102T020000", 15.0, 4)
    prepare(log_path)
    obs = yaml.safe_load(log_path.read_text())["nights"]["2026-01-01"]["observations"]
    assert list(obs) == ["HD-1", "dark", "HD-2", "HD-3"]
    assert obs["HD-3"]["type"] == "star" and obs["HD-3"]["ra_deg"] is None   # unknown to the fake SIMBAD


def test_prepare_unknown_night_without_files(night_dir):
    _, log_path = night_dir
    with pytest.raises(FileNotFoundError):
        prepare_night("2030-01-01", log_path=str(log_path), query_stars=fake_query_stars, log=lambda *_: None)
