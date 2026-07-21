# kbo_occultation/coincidence.py
"""
Multi-telescope coincidence matching for occultation candidates.

A real occultation shadow sweeps across the whole array within
(baseline / shadow velocity) -- milliseconds for the ~100 m MAGIC/LST
baselines -- so genuine events must appear (nearly) simultaneously in
every telescope that was observing, while uncorrelated noise fluctuations 
rarely line up. Requiring a candidate in >= min_ntel telescopes within 
a small time window is therefore a powerful false-positive veto, and its 
residual accidental rate can be measured directly from the data with the 
time-shift method (``accidental_coincidence_rate``).
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .config import ArrayConfig
from .kinematics import ShadowKinematics
from .matched_filter import Candidate, MatchedFilterResult


@dataclass
class TelescopeSearchResult:
    """
    Per-telescope output of a blind matched-filter search, the input
    unit for coincidence matching.
    """
    channel: str
    telescope: str
    candidates: List[Candidate]
    sigma: float
    live_time_s: float
    mf_result: Optional[MatchedFilterResult] = None


@dataclass
class CoincidenceEvent:
    """
    A cluster of per-telescope candidates compatible with one shadow
    crossing.

    time_s : SNR-weighted mean of the member candidate times.
    members : channel -> best (highest-SNR) candidate of that channel
        in the cluster.
    combined_snr : quadrature sum sqrt(sum_i snr_i**2) -- optimal for
        independent Gaussian noise streams.
    """
    time_s: float
    members: Dict[str, Candidate]
    n_tel: int
    combined_snr: float
    max_pairwise_dt_s: float
    max_chi2_reduced: float


@dataclass
class AccidentalRateResult:
    """
    Accidental (chance) coincidence background measured with the
    time-shift method.

    counts : coincidences found in each time-shifted realization.
    rate_per_s : mean accidental rate (counts / live time).
    expected_in_live_time : mean accidental count over the realizations,
        i.e. the expected number of chance coincidences in the actual
        (unshifted) observation.
    """
    counts: np.ndarray
    rate_per_s: float
    expected_in_live_time: float
    tolerance_s: float
    n_shifts: int


def coincidence_tolerance_s(array: ArrayConfig, kinematics: ShadowKinematics,
                            template_duration_s: float,
                            template_allowance: float = 0.5) -> float:
    """
    Coincidence time window: the maximum shadow travel time between
    telescopes plus an allowance for the matched-filter timing jitter
    (a fraction of the template duration, which dominates in practice
    -- ~100 m baselines contribute only milliseconds at ~25 km/s).

    If the shadow drift direction is unknown
    (kinematics.direction_enu is None), the full maximum baseline is
    used isotropically, that is, conservatively.
    """
    if kinematics.direction_enu is None:
        baseline_m = array.max_baseline_m()
    else:
        v_hat = np.asarray(kinematics.direction_enu, dtype=float)
        v_hat = v_hat / np.linalg.norm(v_hat)
        chans = array.channels()
        baseline_m = max(
            abs(float(np.dot(array.baseline_enu_m(a, b), v_hat)))
            for i, a in enumerate(chans) for b in chans[i + 1:]
        )
    return baseline_m / kinematics.v_rel_mps + template_allowance * template_duration_s


def match_coincidences(results: Sequence[TelescopeSearchResult], tolerance_s: float,
                       min_ntel: int = 2,
                       max_pairwise_dt_s: Optional[float] = None) -> List[CoincidenceEvent]:
    """
    Cluster candidates from several telescopes in time and keep the
    clusters seen by at least ``min_ntel`` distinct telescopes.

    Algorithm: pool all candidates tagged by channel, sort by time,
    single-linkage cluster with gap < tolerance_s, then require
    >= min_ntel distinct channels per cluster; within a cluster the
    highest-SNR candidate per channel represents that telescope.

    Parameters
    ----------
    tolerance_s : float
        Single-linkage gap: consecutive candidates closer than this join
        the same cluster.
    min_ntel : int
        Distinct telescopes a cluster must contain to count as an event.
    max_pairwise_dt_s : float, optional
        Extra cut on the *total* span of an event's members (the largest
        time difference between any two of its per-telescope
        candidates). Single-linkage chaining lets a cluster span more
        than tolerance_s when several candidates are strung together, so
        this rejects events whose members can't all belong to one shadow
        crossing. A real crossing sweeps the whole array within
        baseline / v_rel plus timing jitter -- pass that value (e.g. the
        same ``coincidence_tolerance_s`` used for tolerance_s) to enforce
        it. Defaults to None (no span cut -- previous behaviour).
    """
    tagged = [
        (cand.time_s, res.channel, cand)
        for res in results for cand in res.candidates
    ]
    if not tagged:
        return []
    tagged.sort(key=lambda t: t[0])

    events: List[CoincidenceEvent] = []
    cluster = [tagged[0]]                    # start the first cluster with the earliest candidate
    for item in tagged[1:]:
        if item[0] - cluster[-1][0] < tolerance_s:   # item[0] is this candidate's time
            cluster.append(item)             # close enough to the previous one → same cluster
        else:
            ev = _cluster_to_event(cluster, min_ntel, max_pairwise_dt_s)   # gap too big → finalize current cluster
            if ev is not None:
                events.append(ev)
            cluster = [item]                 # start a fresh cluster with this candidate
    ev = _cluster_to_event(cluster, min_ntel, max_pairwise_dt_s)  # to catch the last cluster
    if ev is not None:
        events.append(ev)

    events.sort(key=lambda e: e.combined_snr, reverse=True)
    return events


def _cluster_to_event(cluster, min_ntel: int,
                      max_pairwise_dt_s: Optional[float] = None) -> Optional[CoincidenceEvent]:
    members: Dict[str, Candidate] = {}
    # members dictionary is populated
    for _, channel, cand in cluster:
        if channel not in members or cand.snr > members[channel].snr:
            members[channel] = cand

    if len(members) < min_ntel:
        return None

    # Example: 
    # Cluster = [(100.00, "C", cC), (100.02, "A", cA1), (100.03, "A", cA2), (99.99, "B", cB)]
    # After the loop:
    # members = {"C": cC, "A": (cA1 or cA2, whichever SNR is higher), "B": cB}
    # len(members) == 3

    span = max(c.time_s for c in members.values()) - min(c.time_s for c in members.values())
    if max_pairwise_dt_s is not None and span > max_pairwise_dt_s:
        return None

    snrs = np.array([c.snr for c in members.values()])
    times = np.array([c.time_s for c in members.values()])
    return CoincidenceEvent(
        time_s=float(np.sum(snrs * times) / np.sum(snrs)),
        members=members,
        n_tel=len(members),
        combined_snr=float(np.sqrt(np.sum(snrs**2))),
        max_pairwise_dt_s=float(span),
        max_chi2_reduced=float(max(c.chi2_reduced for c in members.values())),
    )


def accidental_coincidence_rate(results: Sequence[TelescopeSearchResult],
                                tolerance_s: float, live_time_s: float,
                                n_shifts: int = 200, min_shift_s: float = 5.0,
                                min_ntel: int = 2, seed: int = 0,
                                t_start_s: Optional[float] = None,
                                max_pairwise_dt_s: Optional[float] = None) -> AccidentalRateResult:
    """
    Measure the accidental coincidence background with the time-shift method: apply
    a different circular time offset to each telescope's candidate list and re-run
    the coincidence match. Any genuine simultaneous event is destroyed by
    shifts >> tolerance_s, so the coincidences surviving in the shifted realizations
    sample the chance-alignment background of the actual candidate rates (the same
    technique used to background-estimate gravitational-wave coincidence searches).

    Parameters
    ----------
    results : per-telescope search results (the real, unshifted ones).
    tolerance_s, min_ntel, max_pairwise_dt_s : same values used for the
        real match (so the background is estimated under the identical
        cuts).
    live_time_s : duration of the observation the candidates came from.
    n_shifts : number of shifted realizations.
    min_shift_s : minimum relative offset, must be >> tolerance_s.
    t_start_s : start of the circular-shift window. Defaults to the
        earliest candidate time.
    """
    if min_shift_s <= 2 * tolerance_s:
        raise ValueError("min_shift_s must be well above the coincidence tolerance")

    all_times = [c.time_s for res in results for c in res.candidates]
    if t_start_s is None:
        t_start_s = min(all_times) if all_times else 0.0

    rng = np.random.default_rng(seed)
    counts = np.zeros(n_shifts, dtype=int)

    for k in range(n_shifts):
        delta = rng.uniform(min_shift_s, live_time_s - min_shift_s)
        shifted = []
        for j, res in enumerate(results):
            # Telescope j gets offset j * delta (telescope 0 unshifted),
            # wrapped circularly inside the observation window.
            offset = j * delta
            cands = [
                Candidate(
                    time_s=t_start_s + (c.time_s - t_start_s + offset) % live_time_s,
                    snr=c.snr,
                    chi2_reduced=c.chi2_reduced,
                )
                for c in res.candidates
            ]
            shifted.append(TelescopeSearchResult(
                channel=res.channel, telescope=res.telescope, candidates=cands,
                sigma=res.sigma, live_time_s=res.live_time_s,
            ))
        counts[k] = len(match_coincidences(shifted, tolerance_s, min_ntel=min_ntel,
                                           max_pairwise_dt_s=max_pairwise_dt_s))

    mean_count = float(counts.mean())
    return AccidentalRateResult(
        counts=counts,
        rate_per_s=mean_count / live_time_s,
        expected_in_live_time=mean_count,
        tolerance_s=tolerance_s,
        n_shifts=n_shifts,
    )
