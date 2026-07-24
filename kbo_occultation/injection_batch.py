# kbo_occultation/injection_batch.py
"""
Monte Carlo injection-recovery test: does a given noise-reduction method
(highpass, one of the DC-based corrections in dc_combine.py, ...) actually
improve recovery of a real occultation signal, across many KBO radii and
many random injection times -- rather than the one-radius/one-t0 result in
examples/injection_snr_comparison_2025_12_16.py, which turned out to be
dominated by the local noise realization at that single instant rather than
each method's general behavior.

Performance note (see the plan this was built from): a full blind
matched-filter search (matched_filter.sliding_matched_filter_snr) over the
whole light curve costs ~93s per call at real data sizes, almost entirely in
its internal rolling-median baseline (scipy.ndimage.median_filter). Running
that in full for every (radius, t0, method) combination is infeasible (would
be many hours for a useful grid). Two things make this tractable:

1. False-alarm behavior (how many spurious candidates a method's search
   turns up on real noise) doesn't depend on any particular injection, so
   it's measured once per (method, radius) on the true *uninjected* light
   curve (`MethodNoiseProfile`), not repeated per trial.
2. Recovery at a *known* t0 doesn't need a blind search at all -- each
   trial's own corrected+injected series is sliced to a small local window
   around its known t0, and scored there directly. To keep this fast without
   giving up on using each trial's own real (injected) local noise, the
   rolling-median baseline/sigma for a given noise-reduction method is
   computed once (from that method's full, *uninjected*, corrected light
   curve) and reused as an approximation across every trial and radius for
   that method -- justified because a single short injected dip (a fraction
   of a second to a couple of seconds, out of a >7-minute recording) shifts
   a median-filtered baseline by a negligible amount. This is a known,
   accepted approximation, not an exact computation -- see the validation
   spot-check in examples/injection_montecarlo_2025_12_16.py, which compares
   this fast path against a full, uninjected-baseline-free recompute for a
   couple of trials.

Injection itself is always processed one trial at a time: inject, correct,
slice, score, record, discard -- never multiple trials' signals combined
into one array.
"""

from dataclasses import dataclass, asdict
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .config import BandpassConfig, GridConfig, KBOConfig, NumericalConfig, StarConfig
from .simulation import compute_lightcurve
from .detectability import spatial_to_time
from .injection import inject_occultation, random_injection_time
from .matched_filter import compute_baseline_and_sigma, find_candidates, sliding_matched_filter_snr
from .photometry import LightCurve


@dataclass
class TrialResult:
    """One (radius, t0, method) injection trial."""
    radius_m: float
    t0_s: float
    method: str
    recovered: bool
    recovered_snr: float    # nan if nothing survives the shape veto near t0
    offset_ms: float        # nan if not recovered
    chi2_reduced: float     # nan if not recovered
    n_local_candidates: int  # pre-veto candidates found in the local window

    def summary(self) -> dict:
        return asdict(self)


@dataclass
class MethodNoiseProfile:
    """
    One (method, radius) row: false-alarm behavior of a full blind search
    on the true *uninjected* light curve -- computed once, not per trial.
    `channel` is set by the array Monte Carlo (one profile per telescope);
    the single-channel run leaves it empty.
    """
    method: str
    radius_m: float
    sigma: float
    n_candidates: int      # pre-veto
    n_false_alarms: int    # post shape-veto survivors (all false, since uninjected)
    channel: str = ""

    def summary(self) -> dict:
        return asdict(self)


@dataclass
class ArrayTrialResult:
    """
    One (radius, t0, method) injection trial across the whole array: the
    same simulated event injected into every telescope's light curve
    (with the correct ms-scale inter-telescope time offsets), recovered
    per channel, then tested for coincidence.
    """
    radius_m: float
    t0_s: float                          # array-reference injection time
    method: str
    snr_per_channel: Dict[str, float]    # nan where not recovered
    recovered_per_channel: Dict[str, bool]
    n_tel_recovered: int
    coincidence_recovered: bool          # >= min_ntel recoveries within tolerance
    combined_snr: float                  # quadrature over members; nan if no coincidence

    def summary(self) -> dict:
        row = {
            "radius_m": self.radius_m,
            "t0_s": self.t0_s,
            "method": self.method,
            "n_tel_recovered": self.n_tel_recovered,
            "coincidence_recovered": self.coincidence_recovered,
            "combined_snr": self.combined_snr,
        }
        for ch, snr in self.snr_per_channel.items():
            row[f"snr_ch{ch}"] = snr
        for ch, rec in self.recovered_per_channel.items():
            row[f"recovered_ch{ch}"] = rec
        return row


def to_records(items) -> List[dict]:
    """Plain list-of-dicts summary, no pandas required."""
    return [item.summary() for item in items]


def to_dataframe(items):
    """Convert TrialResult/MethodNoiseProfile rows to a pandas DataFrame. Requires pandas."""
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "to_dataframe() requires pandas (pip install pandas); use to_records() instead."
        ) from e
    return pd.DataFrame(to_records(items))


def build_search_template(
    lc: LightCurve,
    correct_fn: Callable[[LightCurve], LightCurve],
    x_m: np.ndarray,
    intensity: np.ndarray,
    shadow_velocity_mps: float,
    template_pad_s: float,
    t0_ref: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    What does this noiseless template look like after passing through
    `correct_fn`, so a matched-filter search on a `correct_fn`-corrected
    light curve compares against a shape that's actually been through the
    same correction -- rather than always matching against the raw,
    unfiltered template. Without this, a correction that reshapes
    sub-cutoff structure (highpass in particular) can still produce a huge
    raw correlation SNR at the true event, yet get rejected by
    matched_filter.shape_veto_chi2's fine-structure check, because the veto
    is comparing filtered data against an unfiltered template -- a test
    mismatch, not evidence the event is undetectable.

    The template is embedded in a realistic, full-duration flat (=1.0)
    baseline on `lc`'s own time grid before correction, rather than
    filtering the short, isolated template array directly -- the latter
    starves a low-cutoff filter of the surrounding context it needs (this
    is exactly the mistake in check_highpass_signal_safety.py's earlier,
    disregarded investigation: a cutoff period comparable to or longer than
    the array being filtered gives meaningless results). Injection is
    placed at the midpoint of `lc.time` purely as a reference location, deep
    in the interior and far from edge effects; matched-filter templates
    only need shape, not absolute time (sliding_matched_filter_snr's
    template_time need not be centered on t=0), so the result is returned
    relative to that reference t0.

    A no-op `correct_fn` (e.g. the "raw" method) reproduces the original
    template almost exactly -- this only matters for corrections that
    reshape structure on the template's own timescale.

    Caution -- this is only valid for corrections that are translation-
    invariant deep in the interior (a true LTI filter like highpass:
    highpass_butterworth(constant) == constant exactly, so a flat reference
    stays flat away from the injected pulse regardless of *where* t0_ref
    is). It is NOT valid for a correction whose reshaping depends on real,
    position-specific auxiliary data -- dc_fft_swap substitutes real DC-
    channel content for any of the template's own sub-cutoff structure, so
    a flat reference comes back showing that specific moment's real trend,
    not a clean, reusable "shape after correction". For such a method, call
    this once per trial with `t0_ref` pinned to that trial's own t0 (so it
    picks up the same real DC value the actual injection would have
    experienced), rather than reusing one canonical call -- see
    run_injection_monte_carlo's `per_trial_reshape_methods`.

    Parameters
    ----------
    lc : LightCurve
        Real light curve, used only for its time grid/length and metadata.
    correct_fn : callable
        The noise-reduction method being tested, same as in method_specs.
    x_m, intensity : array_like
        Noiseless simulated diffraction pattern (as passed to
        injection.inject_occultation).
    shadow_velocity_mps : float
        KBO shadow velocity, for placing the injection.
    template_pad_s : float
        Half-width of the window extracted around the reference t0 --
        should be generous relative to the correction's own memory/ringing
        (e.g. several multiples of 1/cutoff for a highpass filter) while
        staying short enough to still act as a compact matched-filter
        template, not the same as the much larger `pad_s` used for the
        baseline-scoring window elsewhere in this module.
    t0_ref : float, optional
        Reference time to place the injection at. Defaults to the midpoint
        of `lc.time` (deep in the interior, far from edge effects) -- an
        arbitrary but safe choice for a translation-invariant correction.
        Pass the actual trial's t0 for a position-dependent correction.

    Returns
    -------
    search_template_time, search_template_intensity : ndarray
        Time-relative-to-t0_ref template, ready to hand to
        matched_filter.sliding_matched_filter_snr as (template_time,
        template_intensity).
    """
    if t0_ref is None:
        t0_ref = float(lc.time[len(lc.time) // 2])
    flat_lc = LightCurve(lc.time, np.ones_like(lc.signal), meta=dict(lc.meta))
    inj = inject_occultation(flat_lc.time, flat_lc.signal, x_m, intensity, shadow_velocity_mps, t0_s=t0_ref)
    injected_flat_lc = LightCurve(inj.time, inj.injected_signal, meta=dict(flat_lc.meta))
    corrected_flat = correct_fn(injected_flat_lc)

    mask = (corrected_flat.time >= t0_ref - template_pad_s) & (corrected_flat.time <= t0_ref + template_pad_s)
    return corrected_flat.time[mask] - t0_ref, corrected_flat.signal[mask]


def run_injection_monte_carlo(
    lc: LightCurve,
    radii_m: Sequence[float],
    n_trials_per_radius: int,
    grid: GridConfig,
    numerics: NumericalConfig,
    star: StarConfig,
    band: BandpassConfig,
    response: Optional[Callable],
    shadow_velocity_mps: float,
    distance_au: float,
    method_specs: Dict[str, Callable[[LightCurve], LightCurve]],
    pad_s: float,
    snr_threshold: float = 5.0,
    max_chi2_reduced: float = 2.0,
    seed: int = 0,
    checkpoint_path: Optional[str] = None,
    template_pad_s: Optional[float] = None,
    reshape_template_methods: Sequence[str] = (),
    per_trial_reshape_methods: Sequence[str] = (),
) -> Tuple[List[TrialResult], List[MethodNoiseProfile]]:
    """
    Parameters
    ----------
    lc : LightCurve
        Real light curve to inject into (raw, uncorrected).
    radii_m : sequence of float
        KBO radii to test. All must share `grid` (checked below): the
        radius-invariant template duration this relies on only holds for a
        fixed spatial grid.
    n_trials_per_radius : int
        Number of random injection times per radius. The same list of t0s
        (drawn once, see below) is reused across every radius and method.
    grid, numerics, star, band, response, distance_au :
        Passed through to simulation.compute_lightcurve for every radius.
    shadow_velocity_mps : float
        KBO shadow velocity, for spatial_to_time / injection placement.
    method_specs : dict
        {method_name: correction_fn(LightCurve) -> LightCurve}. Include an
        identity function for a "raw" (uncorrected) baseline.
    pad_s : float
        Half-width of the local window scored around each trial's t0. Must
        stay well above sliding_matched_filter_snr's own default detrend
        window (20x template duration) so the local slice's own boundary
        doesn't touch the region actually being scored.
    snr_threshold, max_chi2_reduced : see matched_filter.find_candidates.
    seed : int
        Seeds the shared random t0 list (reproducible across runs).
    checkpoint_path : str, optional
        If given, write the trial results so far to this CSV path after
        each radius completes (requires pandas).
    template_pad_s : float, optional
        Half-width of the per-method "search template" (see
        build_search_template) -- what each method's correction does to the
        template's own shape, used for both correlating and shape-veto
        scoring instead of always assuming the raw, unfiltered template
        applies. Defaults to 10x the (shared) template duration.
    reshape_template_methods : sequence of str
        Names (a subset of method_specs's keys) to build a corrected search
        template for via build_search_template (a single canonical version,
        reused across every radius and trial), instead of using the plain,
        uncorrected template. Only valid for a true LTI filter, translation-
        invariant deep in the interior -- highpass_butterworth(constant) ==
        constant exactly, so a flat reference stays flat away from the
        injected pulse regardless of where the reference injection sits.
        dc_detrend/dc_despike's own memory (10-60s) is two orders of
        magnitude longer than one of these templates' duration (~0.3s), so
        they don't meaningfully reshape sub-event structure either way --
        the plain template remains fine for them (see
        per_trial_reshape_methods for dc_fft_swap, which does reshape but
        isn't translation-invariant). Defaults to an empty tuple (no
        reshaping -- the conservative, always-correct choice); pass the
        highpass-like method names explicitly to opt them in.
    per_trial_reshape_methods : sequence of str
        Names (a subset of method_specs's keys) to build a *fresh* search
        template for on every single trial, pinned to that trial's own t0
        via build_search_template(..., t0_ref=t0_s), instead of one reused
        canonical version. Needed for a correction that reshapes the
        template's own timescale (a genuine chi2 mismatch was measured:
        SNR~165 but chi2/dof~18 for dc_fft_swap at R=300m, the same
        mechanism as highpass) but pulls in real, position-dependent
        auxiliary data -- dc_fft_swap substitutes real DC-channel content
        for any of the template's own sub-cutoff structure, so what a flat
        reference comes back showing depends on *when* you place it (that
        night's real trend at that specific moment), not just the
        correction's fixed impulse response. False-alarm profiling
        (MethodNoiseProfile, which has no single t0 to pin to) still uses
        the plain template for these methods regardless -- it's measuring
        general trigger-happiness on real noise, not this event's shape.

    Returns
    -------
    (trial_results, noise_profiles)
    """
    rng = np.random.default_rng(seed)

    # ─── 1. Simulate the template for every radius; duration must be shared ───
    templates = {}
    for radius_m in radii_m:
        kbo = KBOConfig(radius_m=radius_m, distance_au=distance_au)
        x_m, intensity = compute_lightcurve(kbo, star, band, grid, numerics, response=response)
        template_time = spatial_to_time(x_m, shadow_velocity_mps)
        templates[radius_m] = (x_m, intensity, template_time, float(template_time.max() - template_time.min()))

    durations = {radius_m: templates[radius_m][3] for radius_m in radii_m}
    reference_duration = next(iter(durations.values()))
    for radius_m, duration in durations.items():
        if not np.isclose(duration, reference_duration, rtol=1e-9):
            raise ValueError(
                f"Template duration is not radius-invariant (R={radius_m}m gave {duration}s, "
                f"expected {reference_duration}s) -- the shared-grid assumption this batch "
                f"relies on (one detrend window, one pad_s for every radius) doesn't hold. "
                f"Use a larger/fixed grid.x_max_m so every radius's template fits without "
                f"changing the sampled x_m span."
            )
    common_template_duration_s = reference_duration

    if pad_s < 20.0 * common_template_duration_s:
        raise ValueError(
            f"pad_s={pad_s} is not comfortably above sliding_matched_filter_snr's default "
            f"detrend window (20x template duration = {20.0 * common_template_duration_s}s) -- "
            f"increase pad_s or the local slice's own boundary will bias the score near t0."
        )

    if template_pad_s is None:
        template_pad_s = 10.0 * common_template_duration_s

    # ─── 2. Shared random injection times, reused across every radius/method ───
    any_x_m = templates[radii_m[0]][0]
    t0_list = [
        random_injection_time(lc.time, any_x_m, shadow_velocity_mps, margin_s=pad_s, rng=rng)
        for _ in range(n_trials_per_radius)
    ]

    # ─── 3. Per method: correct the full, uninjected series once; cache baseline/sigma; ───
    # ─── build each method's own "search template" (see build_search_template) ───
    method_cache = {}
    for method, correct_fn in method_specs.items():
        corrected_uninjected = correct_fn(lc)
        baseline, sigma = compute_baseline_and_sigma(
            corrected_uninjected.time, corrected_uninjected.signal, common_template_duration_s
        )
        method_cache[method] = (corrected_uninjected, baseline, sigma)

    # per_trial_reshape_methods use the plain template here -- false-alarm profiling
    # has no single t0 to pin a position-dependent reshape to (see docstring).
    search_templates = {}
    for method, correct_fn in method_specs.items():
        for radius_m in radii_m:
            x_m, intensity, template_time, _ = templates[radius_m]
            if method in reshape_template_methods:
                search_templates[(method, radius_m)] = build_search_template(
                    lc, correct_fn, x_m, intensity, shadow_velocity_mps, template_pad_s
                )
            else:
                search_templates[(method, radius_m)] = (template_time, intensity)

    # ─── 4. Per (method, radius): false-alarm profile on the uninjected series ───
    noise_profiles: List[MethodNoiseProfile] = []
    for method, (corrected_uninjected, baseline, sigma) in method_cache.items():
        for radius_m in radii_m:
            search_time, search_intensity = search_templates[(method, radius_m)]
            mf = sliding_matched_filter_snr(
                corrected_uninjected.time, corrected_uninjected.signal, search_time, search_intensity,
                baseline=baseline, sigma=sigma,
            )
            candidates = find_candidates(mf, snr_threshold=snr_threshold)
            survivors = [c for c in candidates if c.chi2_reduced <= max_chi2_reduced]
            noise_profiles.append(MethodNoiseProfile(
                method=method, radius_m=radius_m, sigma=float(np.median(sigma)),
                n_candidates=len(candidates), n_false_alarms=len(survivors),
            ))

    # ─── 5. Per (radius, t0), one trial at a time: inject, correct, slice, score ───
    trial_results: List[TrialResult] = []
    for radius_m in radii_m:
        x_m, intensity, _, _ = templates[radius_m]

        for t0_s in t0_list:
            inj = inject_occultation(lc.time, lc.signal, x_m, intensity, shadow_velocity_mps, t0_s=t0_s)
            injected_lc = LightCurve(inj.time, inj.injected_signal, meta=dict(lc.meta))

            for method, correct_fn in method_specs.items():
                corrected = correct_fn(injected_lc)
                cached_uninjected, cached_baseline, cached_sigma = method_cache[method]

                if not np.array_equal(corrected.time, cached_uninjected.time):
                    raise AssertionError(
                        f"method {method!r}'s corrected time grid changed between the uninjected "
                        f"reference and this trial's injected series -- the cached-baseline reuse "
                        f"assumes the correction's sample grid is injection-independent (true for "
                        f"despike, which flags only from the untouched DC channel) -- investigate "
                        f"before trusting this method's results."
                    )

                mask = (corrected.time >= t0_s - pad_s) & (corrected.time <= t0_s + pad_s)
                local_time = corrected.time[mask]
                local_signal = corrected.signal[mask]
                local_baseline = cached_baseline[mask]
                local_sigma = cached_sigma[mask]

                if method in per_trial_reshape_methods:
                    # Position-dependent reshaping (e.g. dc_fft_swap): rebuild pinned to
                    # this trial's own t0, so it picks up the real auxiliary data (DC
                    # channel) this specific injection actually experienced.
                    search_time, search_intensity = build_search_template(
                        lc, correct_fn, x_m, intensity, shadow_velocity_mps, template_pad_s, t0_ref=t0_s
                    )
                else:
                    search_time, search_intensity = search_templates[(method, radius_m)]

                mf = sliding_matched_filter_snr(
                    local_time, local_signal, search_time, search_intensity,
                    baseline=local_baseline, sigma=local_sigma,
                )
                candidates = find_candidates(mf, snr_threshold=snr_threshold)
                survivors = [c for c in candidates if c.chi2_reduced <= max_chi2_reduced]

                # A survivor only counts as recovery of *this* injection if it's actually
                # close to t0 -- otherwise "nearest" could be some unrelated false-alarm
                # candidate elsewhere in the local window (same convention as
                # tests/injection_recovery_test.py).
                nearest = min(survivors, key=lambda c: abs(c.time_s - t0_s)) if survivors else None
                if nearest is not None and abs(nearest.time_s - t0_s) < common_template_duration_s:
                    trial_results.append(TrialResult(
                        radius_m=radius_m, t0_s=t0_s, method=method, recovered=True,
                        recovered_snr=nearest.snr, offset_ms=(nearest.time_s - t0_s) * 1e3,
                        chi2_reduced=nearest.chi2_reduced, n_local_candidates=len(candidates),
                    ))
                else:
                    trial_results.append(TrialResult(
                        radius_m=radius_m, t0_s=t0_s, method=method, recovered=False,
                        recovered_snr=np.nan, offset_ms=np.nan, chi2_reduced=np.nan,
                        n_local_candidates=len(candidates),
                    ))

        if checkpoint_path is not None:
            to_dataframe(trial_results).to_csv(checkpoint_path, index=False)

    return trial_results, noise_profiles


def run_array_injection_monte_carlo(
    lcs: Dict[str, LightCurve],
    array,
    radii_m: Sequence[float],
    n_trials_per_radius: int,
    grid: GridConfig,
    numerics: NumericalConfig,
    star: StarConfig,
    band: BandpassConfig,
    response: Optional[Callable],
    kinematics,
    distance_au: float,
    method_specs: Dict[str, Dict[str, Callable[[LightCurve], LightCurve]]],
    pad_s: float,
    snr_threshold: float = 5.0,
    max_chi2_reduced: float = 2.0,
    min_ntel: int = 2,
    max_pairwise_dt_s: Optional[float] = None,
    seed: int = 0,
    checkpoint_path: Optional[str] = None,
    template_pad_s: Optional[float] = None,
    reshape_template_methods: Sequence[str] = (),
    per_trial_reshape_methods: Sequence[str] = (),
) -> Tuple[List["ArrayTrialResult"], List[MethodNoiseProfile]]:
    """
    Array version of run_injection_monte_carlo: the same simulated event
    is injected into every telescope's light curve (offset by the
    shadow's ms-scale travel time between telescopes), each channel is
    scored exactly like the single-channel Monte Carlo, and the trial is
    counted as recovered when >= min_ntel channels recover it within the
    coincidence tolerance. The resulting efficiency curve
    (efficiency_curve) is the one the upper-limit calculation needs.

    Parameters (beyond run_injection_monte_carlo's)
    ----------
    lcs : dict
        channel -> raw LightCurve, e.g. LightCurve.from_stat_binary_all(
        path). Channels not present in `array` are ignored.
    array : config.ArrayConfig
        Telescope positions + channel mapping (e.g. config.ORM_ARRAY).
    kinematics : kinematics.ShadowKinematics or float
        Shadow velocity for the injection stretch and the
        inter-telescope offsets. A bare float means unknown direction
        (offsets = 0, isotropic coincidence tolerance) -- fine, since
        at ~25 km/s the ~100 m baselines contribute only ~4 ms.
    method_specs : dict
        {method: {channel: correction_fn}} -- per-channel correction
        functions, since DC reports differ per telescope. A bare
        callable as a method value is applied to every channel.
    min_ntel : int
        Channels that must recover the injection, within the
        coincidence tolerance, for coincidence_recovered.
    max_pairwise_dt_s : float, optional
        Extra cut on the total span of the recovered channels' times
        (see coincidence.match_coincidences). None uses only the
        single-linkage tolerance.

    Returns
    -------
    (array_trial_results, noise_profiles)
        noise_profiles carry one row per (method, channel, radius).
    """
    from .coincidence import (TelescopeSearchResult, coincidence_tolerance_s,
                              match_coincidences)
    from .kinematics import ShadowKinematics
    from .matched_filter import Candidate

    if not isinstance(kinematics, ShadowKinematics):
        kinematics = ShadowKinematics(v_rel_mps=float(kinematics))
    v_rel = kinematics.v_rel_mps

    channels = [ch for ch in array.channels() if ch in lcs]
    if len(channels) < min_ntel:
        raise ValueError(
            f"Only {len(channels)} of the array's channels are present in lcs "
            f"({channels}), fewer than min_ntel={min_ntel}."
        )

    # Normalize method_specs: allow a single callable per method.
    specs: Dict[str, Dict[str, Callable]] = {}
    for method, spec in method_specs.items():
        if callable(spec):
            specs[method] = {ch: spec for ch in channels}
        else:
            missing = [ch for ch in channels if ch not in spec]
            if missing:
                raise ValueError(f"method {method!r} has no correction for channels {missing}")
            specs[method] = spec

    # Inter-telescope injection offsets: shadow travel time from the
    # reference telescope to each telescope, along the drift direction.
    ref_ch = channels[0]
    if kinematics.direction_enu is None:
        offsets = {ch: 0.0 for ch in channels}
    else:
        v_hat = np.asarray(kinematics.direction_enu, dtype=float)
        v_hat = v_hat / np.linalg.norm(v_hat)
        offsets = {
            ch: float(np.dot(array.baseline_enu_m(ref_ch, ch), v_hat)) / v_rel
            for ch in channels
        }

    rng = np.random.default_rng(seed)

    # ─── 1. Templates per radius; duration must be radius-invariant (shared grid) ───
    templates = {}
    for radius_m in radii_m:
        kbo = KBOConfig(radius_m=radius_m, distance_au=distance_au)
        x_m, intensity = compute_lightcurve(kbo, star, band, grid, numerics, response=response)
        template_time = spatial_to_time(x_m, v_rel)
        templates[radius_m] = (x_m, intensity, template_time, float(template_time.max() - template_time.min()))

    reference_duration = templates[radii_m[0]][3]
    for radius_m in radii_m:
        if not np.isclose(templates[radius_m][3], reference_duration, rtol=1e-9):
            raise ValueError(
                f"Template duration is not radius-invariant (R={radius_m}m gave "
                f"{templates[radius_m][3]}s, expected {reference_duration}s) -- use a "
                f"fixed grid.x_max_m large enough for every radius."
            )
    common_template_duration_s = reference_duration

    if pad_s < 20.0 * common_template_duration_s:
        raise ValueError(
            f"pad_s={pad_s} is not comfortably above the default detrend window "
            f"(20x template duration = {20.0 * common_template_duration_s}s)."
        )
    if template_pad_s is None:
        template_pad_s = 10.0 * common_template_duration_s

    tolerance_s = coincidence_tolerance_s(array, kinematics, common_template_duration_s)

    # ─── 2. Shared random reference injection times (drawn on one channel's grid;
    # all channels share the same time grid from from_stat_binary_all) ───
    any_x_m = templates[radii_m[0]][0]
    t0_list = [
        random_injection_time(lcs[ref_ch].time, any_x_m, v_rel, margin_s=pad_s, rng=rng)
        for _ in range(n_trials_per_radius)
    ]

    # ─── 3. Per (method, channel): correct the uninjected series once, cache
    # baseline/sigma; build search templates ───
    method_cache = {}
    for method in specs:
        for ch in channels:
            corrected_uninjected = specs[method][ch](lcs[ch])
            baseline, sigma = compute_baseline_and_sigma(
                corrected_uninjected.time, corrected_uninjected.signal, common_template_duration_s
            )
            method_cache[(method, ch)] = (corrected_uninjected, baseline, sigma)

    search_templates = {}
    for method in specs:
        for ch in channels:
            for radius_m in radii_m:
                x_m, intensity, template_time, _ = templates[radius_m]
                if method in reshape_template_methods:
                    search_templates[(method, ch, radius_m)] = build_search_template(
                        lcs[ch], specs[method][ch], x_m, intensity, v_rel, template_pad_s
                    )
                else:
                    search_templates[(method, ch, radius_m)] = (template_time, intensity)

    # ─── 4. False-alarm profiles per (method, channel, radius) ───
    noise_profiles: List[MethodNoiseProfile] = []
    for (method, ch), (corrected_uninjected, baseline, sigma) in method_cache.items():
        for radius_m in radii_m:
            search_time, search_intensity = search_templates[(method, ch, radius_m)]
            mf = sliding_matched_filter_snr(
                corrected_uninjected.time, corrected_uninjected.signal, search_time, search_intensity,
                baseline=baseline, sigma=sigma,
            )
            candidates = find_candidates(mf, snr_threshold=snr_threshold)
            survivors = [c for c in candidates if c.chi2_reduced <= max_chi2_reduced]
            noise_profiles.append(MethodNoiseProfile(
                method=method, radius_m=radius_m, sigma=float(np.median(sigma)),
                n_candidates=len(candidates), n_false_alarms=len(survivors),
                channel=ch,
            ))

    # ─── 5. Per (radius, t0, method): inject into every channel, score locally,
    # then test the recoveries for coincidence ───
    trial_results: List[ArrayTrialResult] = []
    for radius_m in radii_m:
        x_m, intensity, _, _ = templates[radius_m]

        for t0_s in t0_list:
            for method in specs:
                snr_per_channel: Dict[str, float] = {}
                recovered_per_channel: Dict[str, bool] = {}
                recovered_candidates: List[TelescopeSearchResult] = []

                for ch in channels:
                    t0_ch = t0_s + offsets[ch]
                    lc = lcs[ch]
                    inj = inject_occultation(lc.time, lc.signal, x_m, intensity, v_rel, t0_s=t0_ch)
                    injected_lc = LightCurve(inj.time, inj.injected_signal, meta=dict(lc.meta))

                    correct_fn = specs[method][ch]
                    corrected = correct_fn(injected_lc)
                    cached_uninjected, cached_baseline, cached_sigma = method_cache[(method, ch)]

                    if not np.array_equal(corrected.time, cached_uninjected.time):
                        raise AssertionError(
                            f"method {method!r} / channel {ch}: corrected time grid changed "
                            f"between the uninjected reference and this trial -- the "
                            f"cached-baseline reuse assumes an injection-independent grid."
                        )

                    mask = (corrected.time >= t0_ch - pad_s) & (corrected.time <= t0_ch + pad_s)
                    local_time = corrected.time[mask]
                    local_signal = corrected.signal[mask]
                    local_baseline = cached_baseline[mask]
                    local_sigma = cached_sigma[mask]

                    if method in per_trial_reshape_methods:
                        search_time, search_intensity = build_search_template(
                            lc, correct_fn, x_m, intensity, v_rel, template_pad_s, t0_ref=t0_ch
                        )
                    else:
                        search_time, search_intensity = search_templates[(method, ch, radius_m)]

                    mf = sliding_matched_filter_snr(
                        local_time, local_signal, search_time, search_intensity,
                        baseline=local_baseline, sigma=local_sigma,
                    )
                    candidates = find_candidates(mf, snr_threshold=snr_threshold)
                    survivors = [c for c in candidates if c.chi2_reduced <= max_chi2_reduced]

                    nearest = min(survivors, key=lambda c: abs(c.time_s - t0_ch)) if survivors else None
                    if nearest is not None and abs(nearest.time_s - t0_ch) < common_template_duration_s:
                        snr_per_channel[ch] = nearest.snr
                        recovered_per_channel[ch] = True
                        recovered_candidates.append(TelescopeSearchResult(
                            channel=ch, telescope=array.by_channel(ch).name,
                            candidates=[nearest], sigma=float(np.median(local_sigma)),
                            live_time_s=float(lc.time[-1] - lc.time[0]),
                        ))
                    else:
                        snr_per_channel[ch] = np.nan
                        recovered_per_channel[ch] = False

                events = match_coincidences(recovered_candidates, tolerance_s, min_ntel=min_ntel,
                                            max_pairwise_dt_s=max_pairwise_dt_s)
                coincidence_recovered = len(events) > 0
                combined_snr = events[0].combined_snr if coincidence_recovered else np.nan

                trial_results.append(ArrayTrialResult(
                    radius_m=radius_m, t0_s=t0_s, method=method,
                    snr_per_channel=snr_per_channel,
                    recovered_per_channel=recovered_per_channel,
                    n_tel_recovered=sum(recovered_per_channel.values()),
                    coincidence_recovered=coincidence_recovered,
                    combined_snr=combined_snr,
                ))

        if checkpoint_path is not None:
            to_dataframe(trial_results).to_csv(checkpoint_path, index=False)

    return trial_results, noise_profiles


def efficiency_curve(trials: Sequence["ArrayTrialResult"], method: str, min_ntel: int = 2):
    """
    Coincidence-level detection efficiency eps(r) for one method, from
    run_array_injection_monte_carlo trials: per radius, the fraction of
    trials recovered in coincidence (as scored by the Monte Carlo's own
    min_ntel). Pass min_ntel=1 to instead count trials recovered by any
    single telescope; other values of min_ntel cannot be re-derived here
    -- rerun the Monte Carlo with the min_ntel you want.

    Returns an upper_limits.EfficiencyCurve, ready to plug into a
    SurveyExposure.
    """
    from .upper_limits import EfficiencyCurve

    rows = [t for t in trials if t.method == method]
    if not rows:
        raise ValueError(f"No trials for method {method!r}")

    radii = np.array(sorted({t.radius_m for t in rows}))
    eff = np.array([
        np.mean([
            (t.coincidence_recovered if min_ntel > 1 else t.n_tel_recovered >= 1)
            for t in rows if t.radius_m == r
        ])
        for r in radii
    ])
    return EfficiencyCurve(radii_m=radii, efficiency=eff)
