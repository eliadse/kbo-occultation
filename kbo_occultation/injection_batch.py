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
    """
    method: str
    radius_m: float
    sigma: float
    n_candidates: int      # pre-veto
    n_false_alarms: int    # post shape-veto survivors (all false, since uninjected)

    def summary(self) -> dict:
        return asdict(self)


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
                method=method, radius_m=radius_m, sigma=sigma,
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
                    baseline=local_baseline, sigma=cached_sigma,
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
