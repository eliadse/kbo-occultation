import pytest

from kbo_occultation.observations.log import load_observation_log, nights_table

LOG = """
site: {extinction_k: {B: 0.23, V: 0.13}}
nights:
  "2026-01-01":
    date: "2026-01-01"
    channels: [A]
    observations:
      HD-1: {file: hd1.npz, type: star}
      HD-2: {file: hd2.npz, type: star, exclude: true}
  "2026-01-02":
    date: "2026-01-02"
    exclude: "no dark run"
    channels: [A]
    observations:
      HD-3: {file: hd3.npz, type: star}
"""


@pytest.fixture
def log(tmp_path, monkeypatch):
    monkeypatch.delenv("KBO_OBSERVATIONS_DIR", raising=False)
    (tmp_path / "observation_log.yaml").write_text(LOG)
    return load_observation_log(str(tmp_path / "observation_log.yaml"))


def test_default_selection_skips_excluded(log):
    nights = log.select()
    assert [n.name for n in nights] == ["2026-01-01"]
    assert list(nights[0].stars) == ["HD-1"]
    assert log.excluded() == ["2026-01-01/HD-2: excluded", "2026-01-02: no dark run"]


def test_named_night_is_used_but_its_excluded_runs_are_not(log):
    assert [n.name for n in log.select(["2026-01-02"])] == ["2026-01-02"]
    assert list(log.select(["2026-01-01"])[0].stars) == ["HD-1"]
    assert log.excluded(["2026-01-01", "2026-01-02"]) == ["2026-01-01/HD-2: excluded"]


def test_include_excluded_keeps_everything(log):
    nights = log.select(include_excluded=True)
    assert [n.name for n in nights] == ["2026-01-01", "2026-01-02"]
    assert list(nights[0].stars) == ["HD-1", "HD-2"]
    # the log itself is never filtered
    assert list(log.nights["2026-01-01"].stars) == ["HD-1", "HD-2"]
    assert nights_table(log)["exclude"].tolist() == [None, "no dark run"]
