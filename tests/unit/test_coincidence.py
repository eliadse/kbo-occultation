import numpy as np
import pytest

from kbo_occultation.config import ORM_ARRAY
from kbo_occultation.analysis.coincidence import (TelescopeSearchResult, accidental_coincidence_rate,
                                         coincidence_tolerance_s, match_coincidences)
from kbo_occultation.simulation.kinematics import ShadowKinematics
from kbo_occultation.analysis.matched_filter import Candidate


def result(channel, times, snrs=None, chi2=1.0, live=450.0):
    snrs = snrs if snrs is not None else [8.0] * len(times)
    cands = [Candidate(time_s=t, snr=s, chi2_reduced=chi2) for t, s in zip(times, snrs)]
    return TelescopeSearchResult(channel=channel, telescope=channel, candidates=cands,
                                 sigma=0.01, live_time_s=live)


def test_threefold_match_within_tolerance():
    res = [
        result("C", [100.0]),
        result("A", [100.02]),
        result("B", [99.99]),
    ]
    events = match_coincidences(res, tolerance_s=0.1, min_ntel=3)
    assert len(events) == 1
    ev = events[0]
    assert ev.n_tel == 3
    assert set(ev.members) == {"A", "B", "C"}
    assert np.isclose(ev.combined_snr, np.sqrt(3 * 8.0**2))
    assert ev.max_pairwise_dt_s <= 0.05


def test_no_match_outside_tolerance():
    res = [result("C", [100.0]), result("A", [101.0]), result("B", [102.0])]
    assert match_coincidences(res, tolerance_s=0.1, min_ntel=2) == []


def test_min_ntel_filtering():
    res = [result("C", [100.0]), result("A", [100.01]), result("B", [200.0])]
    assert len(match_coincidences(res, tolerance_s=0.1, min_ntel=2)) == 1
    assert match_coincidences(res, tolerance_s=0.1, min_ntel=3) == []


def test_best_candidate_per_channel_kept():
    res = [
        result("C", [100.0, 100.03], snrs=[6.0, 9.0]),
        result("A", [100.01], snrs=[7.0]),
    ]
    events = match_coincidences(res, tolerance_s=0.1, min_ntel=2)
    assert len(events) == 1
    assert events[0].members["C"].snr == 9.0


def test_max_pairwise_dt_cut_rejects_chained_cluster():
    # Three candidates each within tolerance of the next, but the two
    # ends span 0.16s -- single-linkage chains them into one cluster.
    res = [
        result("C", [100.00]),
        result("A", [100.08]),
        result("B", [100.16]),
    ]
    # Without the span cut, single-linkage (gap < 0.1) joins all three.
    events = match_coincidences(res, tolerance_s=0.1, min_ntel=3)
    assert len(events) == 1
    assert np.isclose(events[0].max_pairwise_dt_s, 0.16)

    # With a 0.1s span cut, the chained cluster is rejected.
    assert match_coincidences(res, tolerance_s=0.1, min_ntel=3, max_pairwise_dt_s=0.1) == []

    # A genuinely tight coincidence still passes the same cut.
    tight = [result("C", [100.0]), result("A", [100.02]), result("B", [99.99])]
    assert len(match_coincidences(tight, tolerance_s=0.1, min_ntel=3, max_pairwise_dt_s=0.1)) == 1


def test_tolerance_isotropic_and_projected():
    kin_iso = ShadowKinematics(v_rel_mps=25_000.0)
    tol = coincidence_tolerance_s(ORM_ARRAY, kin_iso, template_duration_s=0.3,
                                  template_allowance=0.5)
    expected = ORM_ARRAY.max_baseline_m() / 25_000.0 + 0.15
    assert np.isclose(tol, expected)

    # A drift direction pointing straight up has no baseline projection:
    kin_up = ShadowKinematics(v_rel_mps=25_000.0, direction_enu=np.array([0.0, 0.0, 1.0]))
    tol_up = coincidence_tolerance_s(ORM_ARRAY, kin_up, 0.3, template_allowance=0.5)
    assert np.isclose(tol_up, 0.15)


def test_accidental_rate_matches_analytic_two_fold():
    # Poisson-random candidate lists in two channels: the accidental
    # two-fold rate is R = 2 * r1 * r2 * tau * T events per observation.
    rng = np.random.default_rng(42)
    T = 10_000.0
    n1, n2 = 300, 250
    res = [
        result("C", np.sort(rng.uniform(0, T, n1)), live=T),
        result("A", np.sort(rng.uniform(0, T, n2)), live=T),
    ]
    tau = 0.5
    expected = 2 * (n1 / T) * (n2 / T) * tau * T  # = 7.5

    acc = accidental_coincidence_rate(res, tolerance_s=tau, live_time_s=T,
                                      n_shifts=200, min_shift_s=10.0, seed=1)
    # Statistical agreement within ~4 sigma of the shift-to-shift scatter
    scatter = acc.counts.std() / np.sqrt(acc.n_shifts)
    assert abs(acc.expected_in_live_time - expected) < max(4 * scatter, 0.15 * expected)


def test_accidental_rate_rejects_small_shift():
    res = [result("C", [1.0]), result("A", [2.0])]
    with pytest.raises(ValueError):
        accidental_coincidence_rate(res, tolerance_s=10.0, live_time_s=100.0, min_shift_s=5.0)
