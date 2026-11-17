# kbo_occultation/analysis/noise_removal.py
"""
Which noise-removal method loses the least detection power?

Figure of merit: the detection efficiency eps(R, mag) at a fixed real-data
false-alarm rate. For every method:

1. correct the real run, blind-search it (the same search and empirical
   normalisation as analysis.mimic) and read off the SNR threshold whose
   *real-data* single-telescope FAR equals ``far_per_hr``;
2. inject the star-only occultation (no dimming of the sky background and
   pedestal) at random times into the raw run, apply the method, search with
   the template as the method reshapes it, and count the injections
   recovered above that threshold.

Methods (``build_methods``): Butterworth cutoff scan, PSD-whitened matched
filter, DC detrend window scan, DC FFT swap (crossover <= 0.8 x DC Nyquist,
raised-cosine transition, tapered ends), and the DC despike used as a veto
rather than a correction. Cross-telescope common-mode subtraction is not
offered on purpose: a real event is itself common-mode within <= 8 ms.

Speed trick, exact for every method here: each correction is affine in the
signal (Butterworth, fixed-response whitening, division by a fixed DC trend,
the DC swap). The injected run is therefore
``C(real) + [C(seg * template) - C(seg)]`` on a +-``segment_s`` window, and
only that window needs to be re-corrected and re-searched.
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..config import DC_CUT_HZ, ORM_ARRAY, ArrayConfig
from ..observations.dc_combine import (combine_fast_and_dc, combine_lightcurve_with_dc, detrend_lightcurve_with_dc,
                                       flag_dc_excursions, interpolate_dc_to_fast, load_dc_series,
                                       mask_near_times)
from ..observations.filtering import highpass_lightcurve
from ..observations.photometry import LightCurve, good_time_lightcurve
from ..simulation.kinematics import EARTH_SEMIMAJOR_AXIS_AU, shadow_velocity_opposition
from .injection import inject_occultation
from .mimic import (MimicConfig, _template_duration, bin_lightcurve, search_lightcurve,
                    search_templates, star_bank, star_mas_for)
from .spectra import lightcurve_psd

BUTTERWORTH_CUTOFFS_HZ = (0.3, 0.5, 1.0, 2.0, 3.0)
DC_DETREND_WINDOWS_S = (2.0, 5.0, 10.0)
SWAP_DETREND_WINDOWS_S = (2.0, 3.0, 5.0, 10.0)   # swap first, then DC detrend


@dataclass
class Method:
    """A named correction. ``veto_times`` (unix s) reject candidates within ``veto_pad_s``."""
    name: str
    correct: Callable[[LightCurve], LightCurve]
    veto_times: Optional[np.ndarray] = None
    veto_pad_s: float = 1.0
    needs_dc: bool = False


# ───────────────────────────── whitening ────────────────────────────────────

def whitening_response(lc: LightCurve, n_fft: int, nperseg_target: int = 2 ** 16,
                       smooth_bins: int = 5) -> np.ndarray:
    """
    Gain H(f) (length n_fft//2+1) of the filter that whitens ``lc``'s noise:
    H = sqrt(S_white / S(f)), S_white the median PSD above 100 Hz (or above
    the upper third of the band if the Nyquist is lower), H <= 1 so that the
    white floor passes unchanged and the red noise is suppressed. H(0) = 0
    (the mean is handled by the caller).
    """
    f, psd = lightcurve_psd(lc, nperseg_target)
    psd = np.convolve(psd, np.ones(smooth_bins) / smooth_bins, mode="same") if smooth_bins > 1 else psd
    fs = 1.0 / lc.dt
    floor_sel = f >= min(100.0, 0.66 * fs / 2.0)
    s_white = float(np.median(psd[floor_sel])) if floor_sel.any() else float(np.min(psd[psd > 0]))
    grid = np.fft.rfftfreq(n_fft, d=lc.dt)
    s_grid = np.interp(grid, f, np.maximum(psd, 1e-300))
    h = np.sqrt(np.minimum(1.0, s_white / s_grid))
    h[0] = 0.0
    return h


def make_whitener(reference: LightCurve) -> Callable[[LightCurve], LightCurve]:
    """
    LightCurve -> LightCurve applying the whitening response measured on
    ``reference`` (the real run) to any series on the same cadence. The mean
    is kept, so the result is still a flux near the original level.
    """
    cache: Dict[int, np.ndarray] = {}

    def apply(lc: LightCurve) -> LightCurve:
        n = len(lc.signal)
        if n not in cache:
            cache[n] = whitening_response(reference, n)
        x = np.asarray(lc.signal, dtype=float)
        mean = x.mean()
        y = np.fft.irfft(np.fft.rfft(x - mean) * cache[n], n=n) + mean
        return LightCurve(lc.time, y, meta=dict(lc.meta, correction="psd_whitened"))

    return apply


def swap_then_detrend(lc: LightCurve, dc_report: str, dc_key: str, window_s: float,
                      cut_hz: Optional[float] = None) -> LightCurve:
    """
    DC FFT swap, then division by the smoothed DC trend. ``cut_hz=None`` puts the
    crossover at the Nyquist frequency of this DC series (0.5 / median cadence), bypassing the
    library's 0.8 x Nyquist check on purpose; otherwise it is the given crossover.
    """
    dc_time, dc_value = load_dc_series(dc_report, telescope=dc_key)
    cut = 0.5 / float(np.median(np.diff(dc_time))) if cut_hz is None else cut_hz
    combined = combine_fast_and_dc(lc.time, lc.signal, interpolate_dc_to_fast(lc.time, dc_time, dc_value), cut)
    swapped = LightCurve(lc.time, combined, meta=dict(lc.meta))
    return detrend_lightcurve_with_dc(swapped, dc_report, telescope=dc_key, smooth_window_s=window_s)


# ───────────────────────────── method registry ──────────────────────────────

def build_methods(reference: LightCurve, dc_report: Optional[str], dc_key: str,
                  include: Optional[Sequence[str]] = None) -> Dict[str, Method]:
    """
    All comparison methods for one telescope's run. ``include`` filters by
    name. DC methods are left out when there is no DC report.
    """
    methods: Dict[str, Method] = {"raw": Method("raw", lambda lc: lc)}
    for fc in BUTTERWORTH_CUTOFFS_HZ:
        methods[f"butter_{fc:g}Hz"] = Method(
            f"butter_{fc:g}Hz", lambda lc, fc=fc: highpass_lightcurve(lc, fc))
    methods["whitened"] = Method("whitened", make_whitener(reference))
    if dc_report:
        for w in DC_DETREND_WINDOWS_S:
            methods[f"dc_detrend_{w:g}s"] = Method(
                f"dc_detrend_{w:g}s",
                lambda lc, w=w: detrend_lightcurve_with_dc(lc, dc_report, telescope=dc_key, smooth_window_s=w),
                needs_dc=True)
        for w in SWAP_DETREND_WINDOWS_S:
            methods[f"swapnyq_detrend_{w:g}s"] = Method(
                f"swapnyq_detrend_{w:g}s",
                lambda lc, w=w: swap_then_detrend(lc, dc_report, dc_key, w), needs_dc=True)
        for w in (2.0, 5.0):
            methods[f"swapcfg_detrend_{w:g}s"] = Method(
                f"swapcfg_detrend_{w:g}s",
                lambda lc, w=w: swap_then_detrend(lc, dc_report, dc_key, w, DC_CUT_HZ[dc_key]), needs_dc=True)
        methods["dc_swap"] = Method(
            "dc_swap", lambda lc: combine_lightcurve_with_dc(lc, dc_report, DC_CUT_HZ[dc_key], telescope=dc_key),
            needs_dc=True)
        try:
            t, v = load_dc_series(dc_report, telescope=dc_key, pixel=0)
            sel = (t >= reference.time[0] - 120) & (t <= reference.time[-1] + 120)
            t, v = t[sel], v[sel]
            flagged = t[flag_dc_excursions(t, v)] if len(t) >= 5 else np.array([])
            methods["dc_despike_veto"] = Method("dc_despike_veto", lambda lc: lc, veto_times=flagged,
                                                needs_dc=True)
        except (ValueError, KeyError, IndexError):
            pass
    if include is not None:
        methods = {k: m for k, m in methods.items() if k in include}
    return methods


# ───────────────────────────── core computations ────────────────────────────

@dataclass
class ComparisonConfig:
    mimic: MimicConfig = field(default_factory=MimicConfig)
    far_per_hr: float = 1.0            # real-data single-telescope false alarms per hour
    n_trials: int = 20                 # injections per (run, channel, radius, method)
    segment_s: float = 60.0            # half-width of the re-corrected / re-searched window
    snr_floor: float = 3.0             # candidates below this are never listed
    edge_margin_s: float = 1.0         # array injections stay this far (plus the pattern) inside the good time
    peak_window_s: float = 0.05        # array injections: window around t0 for the peak SNR
    seed: int = 0
    methods: Optional[Sequence[str]] = None


def _crop(lc: LightCurve, a: int, b: int) -> LightCurve:
    return LightCurve(lc.time[a:b], lc.signal[a:b], meta=dict(lc.meta))


def _vetoed(times: np.ndarray, method: Method) -> np.ndarray:
    if method.veto_times is None or len(method.veto_times) == 0:
        return np.zeros(len(times), dtype=bool)
    return mask_near_times(np.asarray(times), method.veto_times, method.veto_pad_s)


def reshaped_templates(lc: LightCurve, method: Method, templates, v_mps: float, cfg: ComparisonConfig):
    """
    {radius: (t, intensity)} as the method reshapes it: the template is
    injected into a flat record, corrected, and the correction of the flat
    record alone is subtracted (so an additive DC trend does not leak in).
    """
    n_half = int(cfg.segment_s / lc.dt)
    mid = len(lc.time) // 2
    seg = _crop(lc, max(0, mid - n_half), min(len(lc.time), mid + n_half))
    flat = LightCurve(seg.time, np.ones_like(seg.signal, dtype=float), meta=dict(lc.meta))
    base = method.correct(flat).signal
    t0 = float(seg.time[len(seg.time) // 2])
    out = {}
    for radius, (tt, ti) in templates.items():
        inj = inject_occultation(flat.time, flat.signal, tt * v_mps, ti, v_mps, t0_s=t0)
        corr = method.correct(LightCurve(inj.time, inj.injected_signal, meta=dict(lc.meta))).signal
        shape = 1.0 + (corr - base) / np.mean(base)
        half = max(abs(tt.min()), abs(tt.max())) + 1.5   # keep the filter's ringing, stay a compact template
        sel = np.abs(flat.time - t0) <= half
        out[radius] = (flat.time[sel] - t0, shape[sel])
    return out


def real_candidates(lc_full: LightCurve, method: Method, templates_m, dur: float, good_time,
                    cfg: ComparisonConfig):
    """Blind search of the corrected, good-time-cropped, binned real run."""
    corrected = bin_lightcurve(good_time_lightcurve(method.correct(lc_full), good_time), cfg.mimic.bin_samples)
    rows, _, _ = search_lightcurve(corrected, templates_m, cfg.snr_floor, dur, cfg.mimic.sigma_window_s)
    rows = [r for r in rows if not _vetoed([r["time_s"]], method)[0]]
    return rows, len(corrected.time) * corrected.dt


def threshold_for_far(snrs: np.ndarray, hours: float, far_per_hr: float, floor: float) -> float:
    """Smallest SNR (>= floor) whose real-data rate is <= far_per_hr."""
    if hours <= 0:
        return float("nan")
    s = np.sort(np.asarray(snrs, dtype=float))
    allowed = int(np.floor(far_per_hr * hours))
    if len(s) <= allowed:
        return floor
    # the (allowed+1)-th largest candidate must lie below the threshold
    return float(np.nextafter(s[len(s) - allowed - 1], np.inf))


def injection_efficiency(lc_full: LightCurve, method: Method, corrected_full: LightCurve, templates_raw,
                         templates_m, v_mps: float, star_fraction: float, good_time,
                         threshold: float, cfg: ComparisonConfig, rng: np.random.Generator,
                         radii: Sequence[float]) -> Dict[float, float]:
    """eps(R): fraction of star-only injections recovered above ``threshold``."""
    dur = _template_duration(templates_raw)
    tol = max(dur, 2.0 * cfg.mimic.bin_samples * lc_full.dt)
    t_lo, t_hi = lc_full.time[0], lc_full.time[-1]
    if good_time is not None:
        t_lo, t_hi = lc_full.time[0] + good_time[0], lc_full.time[0] + good_time[1]
    lo, hi = t_lo + cfg.segment_s, t_hi - cfg.segment_s
    if hi <= lo:
        return {r: float("nan") for r in radii}
    n_half = int(cfg.segment_s / lc_full.dt)
    t0s = rng.uniform(lo, hi, cfg.n_trials)
    hits = {r: 0 for r in radii}
    for t0 in t0s:
        c = int(round((t0 - lc_full.time[0]) / lc_full.dt))
        a, b = max(0, c - n_half), min(len(lc_full.time), c + n_half)
        seg = _crop(lc_full, a, b)
        clean = method.correct(seg).signal
        base_seg = corrected_full.signal[a:b]
        veto = bool(_vetoed([t0], method)[0])
        for radius in radii:
            tt, ti = templates_raw[radius]
            inj = inject_occultation(seg.time, seg.signal, tt * v_mps, ti, v_mps, t0_s=float(t0),
                                     star_fraction=star_fraction)
            corr = method.correct(LightCurve(seg.time, inj.injected_signal, meta=dict(seg.meta))).signal
            window = LightCurve(seg.time, base_seg + (corr - clean), meta=dict(seg.meta))
            window = bin_lightcurve(window, cfg.mimic.bin_samples)
            rows, _, _ = search_lightcurve(window, {radius: templates_m[radius]}, threshold, dur,
                                           cfg.mimic.sigma_window_s)
            ok = any(abs(r["time_s"] - t0) <= tol and r["snr"] >= threshold for r in rows)
            hits[radius] += int(ok and not veto)
    return {r: hits[r] / len(t0s) for r in radii}


# ───────────────────────────── array injections ─────────────────────────────

def shadow_direction_enu(ra_deg: Optional[float], dec_deg: Optional[float], t_unix: float,
                         distance_au: float = 39.0) -> Optional[np.ndarray]:
    """
    Unit ENU vector of the shadow's drift at the site (Earth-reflex convention, the one used for
    coincidence timing), or None when the star's coordinates are unknown.
    """
    if ra_deg is None or dec_deg is None:
        return None
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from ..observations.cc_reports import site_location
    from ..simulation.kinematics import shadow_velocity_from_target
    kin = shadow_velocity_from_target(SkyCoord(ra_deg * u.deg, dec_deg * u.deg), Time(t_unix, format="unix"),
                                      distance_au, include_kbo_orbit=False, location=site_location())
    return kin.direction_enu


def shadow_offsets_s(direction_enu: Optional[np.ndarray], channels: Sequence[str], v_mps: float,
                     rng: Optional[np.random.Generator] = None, array: ArrayConfig = ORM_ARRAY) -> Dict[str, float]:
    """
    Time at which the shadow reaches each telescope relative to the first of ``channels``:
    baseline . v_hat / v. With an unknown direction (None) a horizontal azimuth is drawn from
    ``rng`` (isotropic), as the coincidence tolerance assumes.
    """
    if direction_enu is None:
        az = (rng or np.random.default_rng()).uniform(0.0, 2.0 * np.pi)
        v_hat = np.array([np.sin(az), np.cos(az), 0.0])
    else:
        v_hat = np.asarray(direction_enu, dtype=float)
        v_hat = v_hat / np.linalg.norm(v_hat)
    ref = channels[0]
    return {ch: float(np.dot(array.baseline_enu_m(ref, ch), v_hat)) / v_mps for ch in channels}


def injection_efficiency_array(lcs_full: Dict[str, LightCurve], methods: Dict[str, Method],
                               patterns: Dict[str, Dict[float, Tuple[np.ndarray, np.ndarray]]],
                               templates_m: Dict[str, dict], dur: float, v_mps: float,
                               star_fractions: Dict[str, float], good_time, cfg: ComparisonConfig,
                               rng: np.random.Generator, radii: Sequence[float], *,
                               direction_enu: Optional[np.ndarray] = None,
                               corrected_full: Optional[Dict[str, LightCurve]] = None,
                               array: ArrayConfig = ORM_ARRAY) -> Tuple[List[dict], List[dict]]:
    """
    Star-only injections of the same event into every channel of one run, searched as the real
    data are (``real_candidates``: corrected, good-time cropped, binned, blind search with the
    reshaped template bank), so any threshold and any k-of-n coincidence can be read off afterwards.

    All arguments but ``dur``, ``v_mps``, ``good_time``, ``cfg``, ``rng``, ``radii`` and
    ``direction_enu`` are dicts keyed by channel. ``patterns[ch][R] = (x_m, I)`` is the injected
    pattern (any instrument response), ``templates_m[ch]`` the bank as the method reshapes it,
    ``dur`` the raw template duration of the real search. The same random event times are used
    for all channels; each channel's pattern is shifted by the shadow travel time
    (``shadow_offsets_s``, a new random azimuth per trial when ``direction_enu`` is None).
    Same affine trick as ``injection_efficiency``: only a +-``segment_s`` window is re-corrected.

    Returns ``(trials, candidates)`` row lists. trials: one row per (trial, channel, radius) with
    t0_s, offset_s, star_fraction, peak_snr / peak_dt_s (largest SNR of the bank within
    ``cfg.peak_window_s`` of the event, above the floor or not). candidates: every blind-search
    candidate within max(dur, 2 bins) of the event (time_s, dt_s, snr, chi2_reduced, template radius).
    """
    channels = list(lcs_full)
    corrected_full = corrected_full or {ch: methods[ch].correct(lcs_full[ch]) for ch in channels}
    t_lo = max(lcs_full[ch].time[0] + (good_time[0] if good_time is not None else 0.0) for ch in channels)
    t_hi = min(min(lcs_full[ch].time[0] + good_time[1], lcs_full[ch].time[-1]) if good_time is not None
               else lcs_full[ch].time[-1] for ch in channels)
    half_pattern = max(float(np.abs(x).max()) for ch in channels for x, _ in patterns[ch].values()) / v_mps
    margin = half_pattern + cfg.edge_margin_s
    if t_hi - t_lo <= 2.0 * margin:
        return [], []
    t0s = rng.uniform(t_lo + margin, t_hi - margin, cfg.n_trials)
    bin_n = cfg.mimic.bin_samples
    tol = max(dur, 2.0 * bin_n * lcs_full[channels[0]].dt)
    trials, cands = [], []
    for k, t0 in enumerate(t0s):
        offsets = shadow_offsets_s(direction_enu, channels, v_mps, rng, array)
        for ch in channels:
            lc, method = lcs_full[ch], methods[ch]
            t0c = float(t0 + offsets[ch])
            n_half = int(cfg.segment_s / lc.dt)
            c = int(round((t0 - lc.time[0]) / lc.dt))
            a, b = max(0, c - n_half), min(len(lc.time), c + n_half)
            seg = _crop(lc, a, b)
            clean = method.correct(seg).signal
            base = np.asarray(corrected_full[ch].signal[a:b], dtype=float)
            keep = (seg.time >= t_lo) & (seg.time <= t_hi)
            for radius in radii:
                x_m, intensity = patterns[ch][radius]
                inj = inject_occultation(seg.time, seg.signal, x_m, intensity, v_mps, t0_s=t0c,
                                         star_fraction=star_fractions[ch])
                corr = method.correct(LightCurve(seg.time, inj.injected_signal, meta=dict(seg.meta))).signal
                window = bin_lightcurve(LightCurve(seg.time[keep], (base + (corr - clean))[keep],
                                                   meta=dict(seg.meta)), bin_n)
                rows, _, calib = search_lightcurve(window, templates_m[ch], cfg.snr_floor, dur,
                                                   cfg.mimic.sigma_window_s,
                                                   peak_near=(t0c, cfg.peak_window_s))
                trials.append({"trial": k, "channel": ch, "radius_m": float(radius), "t0_s": t0c,
                               "offset_s": offsets[ch], "star_fraction": star_fractions[ch],
                               "peak_snr": calib["peak_snr"], "peak_dt_s": calib["peak_time"] - t0c})
                cands.extend({"trial": k, "channel": ch, "radius_m": float(radius), "t0_s": t0c,
                              "time_s": r["time_s"], "dt_s": r["time_s"] - t0c, "snr": r["snr"],
                              "chi2_reduced": r["chi2_reduced"], "template_radius_m": r["radius_m"]}
                             for r in rows if abs(r["time_s"] - t0c) <= tol)
    return trials, cands


def star_fraction_for(obs, telescope_type: str = "MAGIC", bin_s: float = 524.288e-6) -> float:
    """S / (S + B) from the instrument model; 1.0 when the magnitude is unknown."""
    if getattr(obs, "mag_B", None) is None:
        return 1.0
    try:
        import astropy.units as u
        from ..simulation.instruments import Instrument
        b = Instrument(telescope_type).noise_budget(obs.mag_B, bin_s * u.s)
        return float(b["star_photons"] / (b["star_photons"] + b["nsb_photons"]))
    except Exception:
        return 1.0


def run_comparison(log, nights=None, cfg: Optional[ComparisonConfig] = None, *,
                   include_excluded: bool = False, progress: Optional[Callable[[str], None]] = None):
    """
    Efficiency at fixed real-data FAR for every method on the selected runs.

    Returns ``(thresholds, efficiency, candidates)`` DataFrames:
    thresholds: one row per method (pooled over runs and channels);
    efficiency: night, observation, channel, mag_B, method, radius_m, efficiency, threshold;
    candidates: every real candidate above ``snr_floor`` per method (with the angular radius ``star_mas`` of the
    template bank of its star, to tell lists made with another bank).
    """
    cfg = cfg or ComparisonConfig()
    say = progress or (lambda m: None)
    rng = np.random.default_rng(cfg.seed)
    v = shadow_velocity_opposition(cfg.mimic.distance_au + EARTH_SEMIMAJOR_AXIS_AU)
    reference = search_templates(cfg.mimic, v)                  # radii and span are those of every star's bank
    radii = list(reference)
    dur = _template_duration(reference)

    jobs = []     # (night, obs, ch, lc_full, methods, templates_m)
    cand_rows, hours_by_method = [], {}
    snr_by_method: Dict[str, List[float]] = {}
    for night in log.select(nights, include_excluded):
        for obs in night.stars.values():
            if not obs.has_cache:
                continue
            templates_raw = star_bank(cfg.mimic, obs, v)          # the bank of this star's size
            star_mas = star_mas_for(cfg.mimic, obs)
            for ch in night.channels:
                lc_full = obs.lightcurve(ch)
                dc_key = night.telescope(ch).dc_key
                methods = build_methods(lc_full, night.dc_report, dc_key, cfg.methods)
                for name, m in methods.items():
                    say(f"{night.name} {obs.name} {ch} {name}: real search")
                    try:
                        tm = reshaped_templates(lc_full, m, templates_raw, v, cfg)
                        rows, live = real_candidates(lc_full, m, tm, dur, obs.good_time, cfg)
                    except (ValueError, KeyError) as e:
                        say(f"  skipped: {e}")
                        continue
                    hours_by_method[name] = hours_by_method.get(name, 0.0) + live / 3600.0
                    snr_by_method.setdefault(name, []).extend(r["snr"] for r in rows)
                    for r in rows:
                        cand_rows.append({"night": night.name, "observation": obs.name, "channel": ch,
                                          "method": name, "star_mas": star_mas, **r})
                    jobs.append((night, obs, ch, lc_full, m, tm, templates_raw))

    thr_rows, thr = [], {}
    for name, hours in hours_by_method.items():
        s = np.array(snr_by_method.get(name, []))
        thr[name] = threshold_for_far(s, hours, cfg.far_per_hr, cfg.snr_floor)
        thr_rows.append({"method": name, "hours": hours, "n_candidates": len(s),
                         "far_per_hr": cfg.far_per_hr, "threshold": thr[name]})

    eff_rows = []
    for night, obs, ch, lc_full, m, tm, templates_raw in jobs:
        say(f"{night.name} {obs.name} {ch} {m.name}: injections")
        corrected_full = m.correct(lc_full)
        sf = star_fraction_for(obs, night.telescope(ch).name.split("-")[0])
        eff = injection_efficiency(lc_full, m, corrected_full, templates_raw, tm, v, sf, obs.good_time,
                                   thr[m.name], cfg, rng, radii)
        for radius, e in eff.items():
            eff_rows.append({"night": night.name, "observation": obs.name, "channel": ch,
                             "mag_B": obs.mag_B, "method": m.name, "radius_m": radius,
                             "efficiency": e, "threshold": thr[m.name], "star_fraction": sf})
    return pd.DataFrame(thr_rows), pd.DataFrame(eff_rows), pd.DataFrame(cand_rows)


def summarize(efficiency: pd.DataFrame) -> pd.DataFrame:
    """Mean efficiency per (method, radius), averaged over runs and channels."""
    if efficiency.empty:
        return efficiency
    return (efficiency.pivot_table(index="method", columns="radius_m", values="efficiency", aggfunc="mean")
            .assign(mean_all=lambda d: d.mean(axis=1)).sort_values("mean_all", ascending=False))
