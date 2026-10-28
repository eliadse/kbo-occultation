import numpy as np
import pandas as pd
import pytest

from kbo_occultation.analysis.mimic import (MimicResult, bin_lightcurve, exceedances, far_table,
                                            search_lightcurve, tail_table)
from kbo_occultation.analysis.coincidence import (TelescopeSearchResult,
                                                  accidental_coincidence_rate)
from kbo_occultation.analysis.matched_filter import Candidate
from kbo_occultation.observations.photometry import LightCurve

DT = 0.005


def _template():
    t = np.arange(-0.4, 0.4, DT)
    return t, 1.0 - 0.25 * np.exp(-0.5 * (t / 0.08) ** 2) * np.cos(8 * t) ** 2


def test_bin_lightcurve_keeps_gaps_and_means():
    t = np.arange(10000) * 0.0005
    keep = np.ones(10000, bool)
    keep[4000:4300] = False
    lc = LightCurve(t[keep], np.arange(10000.0)[keep])
    b = bin_lightcurve(lc, 10)
    assert len(b.time) <= len(lc.time) // 10
    gaps = np.diff(b.time) > 1.5 * 10 * 0.0005
    assert gaps.sum() == 1
    assert b.signal[0] == pytest.approx(np.arange(10).mean())


def test_search_finds_injected_event_with_depth_and_radius():
    rng = np.random.default_rng(0)
    n = 60000
    t = np.arange(n) * DT
    sig = 1.0 + 0.01 * rng.standard_normal(n)
    tt, ti = _template()
    i0 = 30000
    sig[i0 - len(tt) // 2:i0 - len(tt) // 2 + len(tt)] *= ti
    templates = {100.0: (tt, ti)}
    rows, z, calib = search_lightcurve(LightCurve(t, sig), templates, 4.0, 0.8)
    assert rows, "injected event not found"
    top = max(rows, key=lambda r: r["snr"])
    assert abs(top["time_s"] - t[i0]) < 0.05
    assert top["snr"] > 10 and top["radius_m"] == 100.0
    assert top["depth"] == pytest.approx(0.25, rel=0.25)
    assert 0.7 < calib["snr_mad_sigma"] < 1.4


def test_exceedances_counts_runs_not_samples():
    z = np.zeros(100)
    z[10:15] = -4
    z[50] = -4
    z[60:62] = 4
    assert len(exceedances(z, 3, -1)) == 2
    assert len(exceedances(z, 3, +1)) == 1


def _result(cand_rows, event_rows, live_rows):
    return MimicResult(candidates=pd.DataFrame(cand_rows), events=pd.DataFrame(event_rows),
                       live=pd.DataFrame(live_rows), tails=pd.DataFrame(),
                       reference_flags=pd.DataFrame())


def test_far_table_counts_per_telescope_hour():
    key = {"night": "n", "observation": "o", "correction": "raw"}
    cands = [{**key, "kind": "real", "k": -1, "channel": "A", "snr": s, "time_s": i}
             for i, s in enumerate([3.2, 4.1, 5.5, 6.0])]
    cands += [{**key, "kind": "surrogate", "k": 0, "channel": "A", "snr": 3.5, "time_s": 1.0}]
    live = [{**key, "kind": "real", "live_s": 1800.0, "n_tel": 1},
            {**key, "kind": "surrogate", "live_s": 3600.0, "n_tel": 1}]
    far = far_table(_result(cands, [], live), snr_grid=(3, 5))
    r5 = far[(far.kind == "real") & (far.snr_min == 5)].iloc[0]
    assert r5["count"] == 2 and r5["hours"] == pytest.approx(0.5) and r5["rate_per_hr"] == pytest.approx(4.0)
    s3 = far[(far.kind == "surrogate") & (far.snr_min == 3)].iloc[0]
    assert s3["count"] == 1 and s3["rate_per_hr"] == pytest.approx(1.0)


def test_timeslide_reports_min_snr_of_events():
    rng = np.random.default_rng(1)
    results = [TelescopeSearchResult(channel=ch, telescope=ch, sigma=0.0, live_time_s=600.0,
                                     candidates=[Candidate(float(t), 4.0 + j, 1.0)
                                                 for t in rng.uniform(0, 600, 200)])
               for j, ch in enumerate("AB")]
    acc = accidental_coincidence_rate(results, 0.5, 600.0, n_shifts=20, min_shift_s=5.0, min_ntel=2)
    assert len(acc.event_min_snr) == int(acc.counts.sum()) > 0
    assert acc.event_min_snr.min() >= 4.0
