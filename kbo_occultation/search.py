# kbo_occultation/search.py
"""
End-to-end blind search for KBO occultations in one multi-telescope
observation: load all channels of a stat-binary file, apply a noise
reduction, run the per-telescope matched-filter search, and demand
coincidence between telescopes.

This is the entry point for processing real (on-ecliptic) data when it
arrives; on benchmark stars it measures the false-positive behavior of
the full array pipeline.
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Union

import numpy as np

from .config import (ArrayConfig, BandpassConfig, GridConfig, KBOConfig,
                     NumericalConfig, ORM_ARRAY, StarConfig)
from .coincidence import (AccidentalRateResult, CoincidenceEvent,
                          TelescopeSearchResult, accidental_coincidence_rate,
                          coincidence_tolerance_s, match_coincidences)
from .detectability import spatial_to_time
from .filtering import highpass_lightcurve
from .injection_batch import build_search_template
from .kinematics import ShadowKinematics, shadow_velocity_opposition
from .matched_filter import (Candidate, compute_baseline_and_sigma,
                             find_candidates, sliding_matched_filter_snr)
from .photometry import LightCurve
from .simulation import compute_lightcurve


@dataclass
class ObservationSearchResult:
    per_telescope: Dict[str, TelescopeSearchResult]
    coincidences: List[CoincidenceEvent]
    tolerance_s: float
    accidental: Optional[AccidentalRateResult]
    live_time_s: float
    meta: dict = field(default_factory=dict)


def _load_all_channels(path: str) -> Dict[str, LightCurve]:
    """
    Load every channel of one observation as a dict of LightCurves,
    dispatching on the file type: a pre-processed ``.npz`` cache goes through
    ``from_preprocessed_all`` (raw variance, to match the binary path), anything
    else is parsed as a raw stat-binary via ``from_stat_binary_all``.
    """
    if str(path).endswith(".npz"):
        return LightCurve.from_preprocessed_all(path)
    return LightCurve.from_stat_binary_all(path)


def _dedupe_candidates(candidates: List[Candidate], min_separation_s: float) -> List[Candidate]:
    """
    Merge candidates found by different radius templates at (nearly) the
    same time: keep the highest-SNR one within each min_separation_s
    neighborhood.
    """
    kept: List[Candidate] = []
    for cand in sorted(candidates, key=lambda c: c.snr, reverse=True):
        if all(abs(cand.time_s - k.time_s) >= min_separation_s for k in kept):
            kept.append(cand)
    return kept


def _build_correction(correction: str, dc_reports: Optional[Dict[str, str]],
                      dc_key: str, highpass_cutoff_hz: float) -> Callable[[LightCurve], LightCurve]:
    """Resolve a named correction to a LightCurve -> LightCurve callable."""
    if correction == "raw":
        return lambda lc: lc
    if correction == "highpass":
        return lambda lc: highpass_lightcurve(lc, highpass_cutoff_hz)

    if correction in ("dc_detrend", "dc_despike"):
        if dc_reports is None or dc_key not in dc_reports:
            raise ValueError(
                f"correction {correction!r} needs dc_reports[{dc_key!r}] (path to that "
                f"telescope's DC report pickle)."
            )
        from .dc_combine import despike_lightcurve_with_dc, detrend_lightcurve_with_dc
        path = dc_reports[dc_key]
        if correction == "dc_detrend":
            return lambda lc: detrend_lightcurve_with_dc(lc, path, telescope=dc_key)
        return lambda lc: despike_lightcurve_with_dc(lc, path, telescope=dc_key)[0]

    raise ValueError(
        f"Unknown correction {correction!r}; use 'raw', 'highpass', 'dc_detrend', "
        f"'dc_despike', or pass a callable/dict."
    )


def search_observation(
    bin_file: str,
    dc_reports: Optional[Dict[str, str]] = None,
    array: ArrayConfig = ORM_ARRAY,
    template_radii_m: Sequence[float] = (300.0, 500.0, 1000.0),
    distance_au: float = 43.0,
    star: Optional[StarConfig] = None,
    band: Optional[BandpassConfig] = None,
    grid: Optional[GridConfig] = None,
    numerics: Optional[NumericalConfig] = None,
    response: Optional[Callable] = None,
    kinematics: Optional[Union[ShadowKinematics, float]] = None,
    correction: Union[str, Callable, Dict[str, Callable]] = "raw",
    highpass_cutoff_hz: float = 1.0,
    snr_threshold: float = 5.0,
    max_chi2_reduced: float = 2.0,
    min_ntel: int = 2,
    max_pairwise_dt_s: Optional[float] = None,
    n_background_shifts: int = 200,
) -> ObservationSearchResult:
    """
    Blind coincidence search over one observation.

    Flow: load all channels -> per telescope, apply `correction` ->
    per template radius, sliding matched filter + chi2 shape veto ->
    merge candidates across radii -> demand >= min_ntel telescopes
    within the coincidence tolerance -> measure the accidental
    coincidence background with time shifts.

    Parameters
    ----------
    bin_file : path to the stat-binary observation file (channels A/B/C).
    dc_reports : dict, optional
        dc_key ("m1"/"m2"/"lst1") -> path to that telescope's DC report
        pickle. Only needed for the DC-based corrections.
    array : telescope positions + channel mapping (default ORM_ARRAY).
    template_radii_m : KBO radii to build search templates for.
    star, band, grid, numerics, response :
        Simulation configuration for the templates. Defaults: a
        sun-like point-ish star (0.01 mas), 300-650 nm, a 20 km grid.
    kinematics : ShadowKinematics or float, optional
        Shadow velocity (and optionally drift direction) for template
        stretching and the coincidence tolerance. Default: opposition
        velocity at `distance_au`, unknown direction (isotropic,
        conservative tolerance).
    correction : str, callable, or {channel: callable}
        "raw", "highpass", "dc_detrend", "dc_despike", a single
        LightCurve -> LightCurve callable applied to every channel, or
        a per-channel dict of callables.
    min_ntel : telescopes required in coincidence (default 2 of 3).
    max_pairwise_dt_s : float, optional
        Reject coincidences whose members span more than this (see
        coincidence.match_coincidences). None uses only the
        single-linkage tolerance. Pass ``tolerance_s`` (available on the
        returned result) or a fraction of it to require all members
        within one shadow-crossing window.
    n_background_shifts : time-shift realizations for the accidental
        rate (0 disables the background estimate).

    Returns
    -------
    ObservationSearchResult
    """
    star = star or StarConfig(temperature_K=6000.0, angular_radius_mas=0.01)
    band = band or BandpassConfig(300.0, 650.0, 25)
    grid = grid or GridConfig(x_max_m=20000.0, n_x=2000)
    numerics = numerics or NumericalConfig()

    if kinematics is None:
        kinematics = ShadowKinematics(v_rel_mps=shadow_velocity_opposition(distance_au))
    elif not isinstance(kinematics, ShadowKinematics):
        kinematics = ShadowKinematics(v_rel_mps=float(kinematics))
    v_rel = kinematics.v_rel_mps

    # ─── Load all channels ───
    lcs = _load_all_channels(bin_file)
    channels = [ch for ch in array.channels() if ch in lcs]
    if len(channels) < min_ntel:
        raise ValueError(f"Only channels {channels} available; min_ntel={min_ntel}")

    # ─── Resolve corrections per channel ───
    correct_fns: Dict[str, Callable] = {}
    for ch in channels:
        tel = array.by_channel(ch)
        if callable(correction):
            correct_fns[ch] = correction
        elif isinstance(correction, dict):
            correct_fns[ch] = correction[ch]
        else:
            correct_fns[ch] = _build_correction(correction, dc_reports, tel.dc_key,
                                                highpass_cutoff_hz)

    # ─── Simulate search templates (shared across channels) ───
    templates = {}
    for radius_m in template_radii_m:
        kbo = KBOConfig(radius_m=radius_m, distance_au=distance_au)
        x_m, intensity = compute_lightcurve(kbo, star, band, grid, numerics, response=response)
        template_time = spatial_to_time(x_m, v_rel)
        templates[radius_m] = (x_m, intensity, template_time)
    template_duration_s = float(
        templates[template_radii_m[0]][2].max() - templates[template_radii_m[0]][2].min()
    )

    tolerance_s = coincidence_tolerance_s(array, kinematics, template_duration_s)

    # ─── Per-telescope blind search ───
    per_telescope: Dict[str, TelescopeSearchResult] = {}
    for ch in channels:
        tel = array.by_channel(ch)
        corrected = correct_fns[ch](lcs[ch])

        baseline, sigma = compute_baseline_and_sigma(
            corrected.time, corrected.signal, template_duration_s
        )

        all_candidates: List[Candidate] = []
        for radius_m in template_radii_m:
            x_m, intensity, template_time = templates[radius_m]
            if correction == "highpass":
                # LTI correction reshapes the template's own fine
                # structure; match against the corrected shape so the
                # chi2 veto compares like with like.
                search_time, search_intensity = build_search_template(
                    lcs[ch], correct_fns[ch], x_m, intensity, v_rel,
                    template_pad_s=10.0 * template_duration_s,
                )
            else:
                search_time, search_intensity = template_time, intensity

            mf = sliding_matched_filter_snr(
                corrected.time, corrected.signal, search_time, search_intensity,
                baseline=baseline, sigma=sigma,
            )
            all_candidates.extend(find_candidates(
                mf, snr_threshold=snr_threshold, max_chi2_reduced=max_chi2_reduced
            ))

        candidates = _dedupe_candidates(all_candidates, template_duration_s)
        per_telescope[ch] = TelescopeSearchResult(
            channel=ch, telescope=tel.name, candidates=candidates, sigma=sigma,
            live_time_s=float(corrected.time[-1] - corrected.time[0]),
        )

    live_time_s = max(res.live_time_s for res in per_telescope.values())

    # ─── Coincidence + accidental background ───
    # TODO(decide): should max_pairwise_dt_s default to tolerance_s here?
    # tolerance_s is the single-linkage gap, so a chain of candidates can
    # currently span up to (n_members - 1) * tolerance_s. Setting the span
    # cut to tolerance_s would force all members inside one shadow-crossing
    # window -- physically stricter, but check on real 3-telescope data
    # whether it costs genuine 3-fold efficiency before making it the
    # default (measure the accidental rate both ways with n_background_shifts).
    coincidences = match_coincidences(list(per_telescope.values()), tolerance_s,
                                      min_ntel=min_ntel, max_pairwise_dt_s=max_pairwise_dt_s)

    accidental = None
    if n_background_shifts > 0:
        min_shift_s = max(5.0, 3.0 * tolerance_s)
        accidental = accidental_coincidence_rate(
            list(per_telescope.values()), tolerance_s, live_time_s,
            n_shifts=n_background_shifts, min_shift_s=min_shift_s, min_ntel=min_ntel,
            max_pairwise_dt_s=max_pairwise_dt_s,
        )

    return ObservationSearchResult(
        per_telescope=per_telescope,
        coincidences=coincidences,
        tolerance_s=tolerance_s,
        accidental=accidental,
        live_time_s=live_time_s,
        meta={
            "bin_file": bin_file,
            "correction": correction if isinstance(correction, str) else "custom",
            "template_radii_m": tuple(template_radii_m),
            "distance_au": distance_au,
            "v_rel_mps": v_rel,
            "snr_threshold": snr_threshold,
            "max_chi2_reduced": max_chi2_reduced,
            "min_ntel": min_ntel,
        },
    )
