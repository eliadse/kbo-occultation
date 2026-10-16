import datetime as dt
import os

import numpy as np
import yaml

from kbo_occultation.observations.log import load_observation_log, validate
from kbo_occultation.observations.prepare import (build_caches, night_files, parse_filename,
                                                  prepare_night, scaffold_night)

from test_photometry_loader import write_bin

STEM = "Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_{}_10002_20002_{}"


def make_night_dir(tmp_path, runs):
    """Write tiny .bin files named like the real ones; runs = [(name, timestamp)]."""
    rng = np.random.default_rng(0)
    for name, ts in runs:
        write_bin(tmp_path / (STEM.format(name, ts) + ".bin"), rng.normal(3.0, 0.05, 500).astype(np.float32))
    return tmp_path


def test_parse_filename():
    info = parse_filename(STEM.format("HD-4335", "20251216T212254") + ".bin")
    assert info["name"] == "HD-4335"
    assert info["start_utc"] == dt.datetime(2025, 12, 16, 21, 22, 54)
    assert info["daq"]["sample_duration_ns"] == 2
    assert info["daq"]["range_mV"] == 200 and info["daq"]["impedance_ohm"] == 50
    assert parse_filename("something_else.bin") is None


def test_night_window_spans_midnight(tmp_path):
    make_night_dir(tmp_path, [("HD1", "20260819T230000"), ("DarkPatch5", "20260820T005147"),
                              ("HD2", "20260820T130000"), ("HD0", "20260819T110000")])
    names = [parse_filename(p)["name"] for p in night_files("2026-08-19", str(tmp_path))]
    assert names == ["HD1", "DarkPatch5"]


def test_scaffold_is_valid_yaml_and_loads(tmp_path):
    make_night_dir(tmp_path, [("HD-1", "20260101T220000"), ("dark", "20260101T230000")])
    os.makedirs(tmp_path / "DCs" / "2026_01_01")
    (tmp_path / "DCs" / "2026_01_01" / "dc_report.pkl").write_bytes(b"")
    text = scaffold_night("2026-01-01", str(tmp_path))

    block = yaml.safe_load("nights:\n" + text)["nights"]["2026-01-01"]
    assert block["dark_run"] == "dark"
    assert block["dc_report"] == "DCs/2026_01_01/dc_report.pkl"
    assert block["observations"]["HD-1"]["type"] == "star"
    assert block["observations"]["dark"]["type"] == "dark"
    assert block["observations"]["HD-1"]["ra_deg"] is None

    # The stub pasted into a log is loadable as-is and validate() lists the TODOs.
    log_path = tmp_path / "observation_log.yaml"
    log_path.write_text("site: {extinction_k: {B: 0.23, V: 0.13}}\nnights:\n" + text)
    log = load_observation_log(str(log_path))
    problems = validate(log.nights["2026-01-01"])
    assert any("HD-1: no .npz cache yet" in p for p in problems)
    assert any("HD-1: null ra_deg" in p for p in problems)


def test_build_caches_is_idempotent(tmp_path, monkeypatch):
    monkeypatch.delenv("KBO_OBSERVATIONS_DIR", raising=False)
    make_night_dir(tmp_path, [("HD-1", "20260101T220000")])
    log_path = tmp_path / "observation_log.yaml"
    log_path.write_text("site: {}\nnights:\n" + scaffold_night("2026-01-01", str(tmp_path)))
    night = load_observation_log(str(log_path)).nights["2026-01-01"]

    assert len(build_caches(night, log=lambda *_: None)) == 1
    assert night.observations["HD-1"].has_cache
    assert build_caches(night, log=lambda *_: None) == []          # kept
    assert len(build_caches(night, force=True, log=lambda *_: None)) == 1

    # no CC reports listed -> conditions skipped, validate() still runs
    problems = prepare_night(night, log=lambda *_: None)
    assert not any("cache" in p for p in problems)
