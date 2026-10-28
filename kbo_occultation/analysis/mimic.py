# kbo_occultation/analysis/mimic.py
"""
How often does the real noise look like an occultation?

For every logged star run (ObservationLog.select), blind-search the light
curves under several corrections and keep every candidate above a low SNR
threshold, with the auxiliary information that could explain it. Repeat on
Gaussian surrogates (phase-randomised copies: same spectrum, same gaps, no
transients) and on time-shifted coincidences, and compare:

- real FAR above surrogate FAR  -> transient, non-Gaussian mimics (gusts,
  glitches, tracking jumps, transparency changes);
- real coincident FAR above the time-slide FAR -> noise common to the
  telescopes (clouds, atmosphere), not chance alignment.

Entry points: ``run_mimic`` (data -> tables), ``far_table`` (tables ->
events per hour above SNR), ``tail_table``, ``flag_enrichment_table``;
``kbo-night mimic`` writes them to ``outputs/mimic/``.

Details that matter when reading the numbers:

- light curves are block-averaged by ``bin_samples`` (default 10 = 5 ms,
  well above the ~30 Hz signal-band Nyquist) before searching;
- the SNR is the empirical matched-filter SNR (matched_filter,
  normalization="empirical"), so red noise is not under-counted;
- the chi2/dof stored per candidate is not used as a veto here (a veto
  would hide exactly the mimics we want to count);
- the time-slide background is only defined for coincidences. A single
  telescope's candidate count does not change under a time shift.
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.special import ndtr

from ..config import (BandpassConfig, GridConfig, KBOConfig, NumericalConfig,
                      ORM_ARRAY, StarConfig)
from ..observations.cc_reports import load_conditions
from ..observations.dc_combine import flag_dc_excursions, load_dc_series
from ..observations.photometry import LightCurve, good_time_lightcurve
from ..simulation.engine import OccultationEngine
from ..simulation.kinematics import (EARTH_SEMIMAJOR_AXIS_AU, ShadowKinematics,
                                     shadow_velocity_opposition)
from ..simulation.physics import AU_m, mas_to_rad
from .coincidence import (TelescopeSearchResult, accidental_coincidence_rate,
                          coincidence_tolerance_s, match_coincidences)
from .matched_filter import (Candidate, compute_baseline_and_sigma, find_candidates,
                             sliding_matched_filter_snr)
from .search import _build_correction
from .surrogates import phase_randomized_lightcurve

CORRECTIONS = ("raw", "dc_detrend", "dc_despike")
DC_CORRECTIONS = ("dc_detrend", "dc_despike")
TAIL_SIGMAS = (3, 5, 7)


@dataclass
class MimicConfig:
    corrections: Sequence[str] = CORRECTIONS
    snr_threshold: float = 3.0
    template_radii_m: Sequence[float] = (50.0, 100.0, 300.0, 1000.0)
    star_mas: float = 0.03
    distance_au: float = 40.0
    star_temperature_K: float = 10000.0
    band: BandpassConfig = field(default_factory=lambda: BandpassConfig(320.0, 450.0, 21))
    bin_samples: int = 10
    n_surrogates: int = 20
    n_shifts: int = 200
    aux_window_s: float = 1.0
    n_reference_times: int = 2000
    seed: int = 0
    sigma_window_s: Optional[float] = None   # window of the empirical filter-output sigma (None = detrend window)
    pixel_target: int = 0
    pixels_background: Tuple[int, ...] = (1, 2)


# ───────────────────────── templates and searching ──────────────────────────

def search_templates(cfg: MimicConfig, v_mps: float) -> Dict[float, Tuple[np.ndarray, np.ndarray]]:
    """{radius: (template_time_s, intensity)}, shared span so one detrend window fits all."""
    half = 10000.0
    grid = GridConfig(half, int(2 * half / 5.0))
    r_star = cfg.star_mas * mas_to_rad * cfg.distance_au * AU_m
    engine = OccultationEngine(StarConfig(cfg.star_temperature_K, cfg.star_mas), cfg.band, grid,
                               NumericalConfig(n_r_grid=4000, n_star_side=48),
                               r_max_m=1.3 * (half + 2 * r_star))
    out = {}
    for radius in cfg.template_radii_m:
        x_m, intensity = engine.compute(KBOConfig(radius, cfg.distance_au, 0.0), cfg.star_mas)
        out[float(radius)] = (x_m / v_mps, intensity)
    return out


def bin_lightcurve(lc: LightCurve, n: int) -> LightCurve:
    """
    Block-average ``n`` samples. The series is split where the time step
    exceeds 1.5 dt (removed samples), each run binned on its own, so bins
    never straddle a gap and the gap mask survives.
    """
    if n <= 1:
        return lc
    time, signal = np.asarray(lc.time), np.asarray(lc.signal)
    breaks = np.where(np.diff(time) > 1.5 * lc.dt)[0] + 1
    bounds = np.concatenate([[0], breaks, [len(time)]])
    ts, ss = [], []
    for a, b in zip(bounds[:-1], bounds[1:]):
        m = (b - a) // n
        if m < 1:
            continue
        ts.append(time[a:a + m * n].reshape(m, n).mean(axis=1))
        ss.append(signal[a:a + m * n].reshape(m, n).mean(axis=1))
    return LightCurve(np.concatenate(ts), np.concatenate(ss), meta=dict(lc.meta, binned=n))


def _template_duration(templates) -> float:
    t = next(iter(templates.values()))[0]
    return float(t.max() - t.min())


def search_lightcurve(lc: LightCurve, templates, snr_threshold: float,
                      template_duration_s: float, sigma_window_s: Optional[float] = None):
    """
    Blind search of one (binned, corrected) light curve with every template.
    Returns (rows, z, calib): candidate dicts (best template per time), the
    normalised residual used for the tail statistics, and the calibration of
    the SNR trace itself (std and MAD-sigma, averaged over templates; both
    should be ~1 on noise the normalisation describes correctly).
    """
    baseline, sigma = compute_baseline_and_sigma(lc.time, lc.signal, template_duration_s)
    best: List[dict] = []
    stds, mads = [], []
    for radius, (tt, ti) in templates.items():
        mf = sliding_matched_filter_snr(lc.time, lc.signal, tt, ti, baseline=baseline, sigma=sigma,
                                        normalization="empirical", sigma_window_s=sigma_window_s)
        stds.append(float(np.std(mf.snr)))
        mads.append(float(1.4826 * np.median(np.abs(mf.snr - np.median(mf.snr)))))
        peak_amp = float(np.max(np.abs(mf.residual_template)))
        for cand in find_candidates(mf, snr_threshold=snr_threshold, chi2_weighting="empirical"):
            lag = int(np.searchsorted(mf.time, cand.time_s))
            lag = min(lag, len(mf.amp_sigma) - 1)
            amp = cand.snr * mf.amp_sigma[lag] / mf.template_norm
            best.append({"time_s": cand.time_s, "snr": cand.snr, "chi2_reduced": cand.chi2_reduced,
                         "radius_m": radius, "depth": float(amp * peak_amp)})
    kept: List[dict] = []
    for row in sorted(best, key=lambda r: -r["snr"]):
        if all(abs(row["time_s"] - k["time_s"]) >= template_duration_s for k in kept):
            kept.append(row)
    z = (lc.signal / baseline - 1.0) / np.where(sigma > 0, sigma, np.nan)
    return sorted(kept, key=lambda r: r["time_s"]), z, {"snr_std": float(np.mean(stds)),
                                                         "snr_mad_sigma": float(np.mean(mads))}


def _as_candidates(rows: Iterable[dict]) -> List[Candidate]:
    return [Candidate(time_s=r["time_s"], snr=r["snr"], chi2_reduced=r["chi2_reduced"]) for r in rows]


def _telescope_results(per_channel: Dict[str, Tuple[List[dict], float]], array=ORM_ARRAY):
    return [TelescopeSearchResult(channel=ch, telescope=array.by_channel(ch).name,
                                  candidates=_as_candidates(rows), sigma=0.0, live_time_s=live)
            for ch, (rows, live) in per_channel.items()]


def _event_rows(results: Sequence[TelescopeSearchResult], tolerance_s: float, min_ntel: int = 2):
    events = match_coincidences(list(results), tolerance_s, min_ntel=min_ntel)
    return [{"time_s": ev.time_s, "n_tel": ev.n_tel, "combined_snr": ev.combined_snr,
             "min_snr": min(c.snr for c in ev.members.values()),
             "max_chi2_reduced": ev.max_chi2_reduced,
             "channels": "".join(sorted(ev.members))} for ev in events]


def exceedances(z: np.ndarray, k: float, sign: int = -1) -> np.ndarray:
    """Indices where each contiguous run of ``sign * z > k`` starts."""
    flag = np.nan_to_num(sign * z, nan=0.0) > k
    return np.flatnonzero(flag & ~np.concatenate([[False], flag[:-1]]))


# ───────────────────────── auxiliary information ────────────────────────────

class AuxData:
    """
    Everything that could explain a candidate at time t on a given channel,
    loaded once per night: DC on-target and background pixels with their
    excursion flags, starguider deviation, wind.

    Cadences: DC 1-2 s, starguider ~10 s, weather ~16 s. A +-1 s window
    contains one or two DC samples (the window is widened to hold >= 3) and
    usually no starguider / weather report, so those two are linearly
    interpolated at t, and their run-level statistics are kept as
    covariates.
    """

    def __init__(self, night, cfg: MimicConfig, t0: float, t1: float):
        self.cfg = cfg
        self.t0, self.t1 = t0, t1
        self.night = night
        self.cond = load_conditions(night.conditions_dir) if night.conditions_dir else {}
        self._dc: Dict[Tuple[str, int], Optional[tuple]] = {}

    def _dc_series(self, dc_key: str, pixel: int):
        key = (dc_key, pixel)
        if key not in self._dc:
            out = None
            if self.night.dc_report:
                try:
                    t, v = load_dc_series(self.night.dc_report, telescope=dc_key, pixel=pixel)
                    pad = 120.0
                    sel = (t >= self.t0 - pad) & (t <= self.t1 + pad)
                    t, v = t[sel], v[sel]
                    if len(t) >= 5:
                        out = (t, v, flag_dc_excursions(t, v), float(np.median(np.diff(t))))
                except (ValueError, KeyError, AttributeError, IndexError):
                    out = None
            self._dc[key] = out
        return self._dc[key]

    def _interp(self, tel_key: str, report: str, columns: Sequence[str], t):
        df = self.cond.get(tel_key, {}).get(report)
        if df is None or df.empty or not set(columns) <= set(df.columns):
            return {c: np.full(np.shape(t), np.nan) for c in columns}
        d = df.sort_values("time")
        return {c: np.interp(t, d["time"].to_numpy(), d[c].to_numpy(), left=np.nan, right=np.nan)
                for c in columns}

    def features(self, channel: str, t) -> Dict[str, np.ndarray]:
        """Auxiliary features at times ``t`` (array, unix s) of ``channel``."""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        dc_key = self.night.telescope(channel).dc_key
        n = len(t)
        out = {k: np.full(n, np.nan) for k in
               ("dc_slope_rel_per_s", "dc_excursion", "bkg_dc_excursion",
                "starg_dev_deg", "wind_kmh", "wind_peak_kmh")}
        win = self.cfg.aux_window_s

        on = self._dc_series(dc_key, self.cfg.pixel_target)
        if on is not None:
            dct, dcv, bad, dt = on
            w = max(win, 1.5 * dt)
            bad_t = dct[bad]
            out["dc_excursion"] = _near(t, bad_t, w + dt / 2)
            slope = np.full(n, np.nan)
            for i, ti in enumerate(t):
                s = (dct >= ti - w) & (dct <= ti + w)
                if s.sum() >= 3:
                    slope[i] = np.polyfit(dct[s] - ti, dcv[s], 1)[0] / max(abs(dcv[s].mean()), 1e-12)
            out["dc_slope_rel_per_s"] = slope
        bkg = [self._dc_series(dc_key, p) for p in self.cfg.pixels_background]
        bkg = [b for b in bkg if b is not None]
        if bkg:
            hit = np.zeros(n, dtype=bool)
            for dct, _, bad, dt in bkg:
                hit |= _near(t, dct[bad], max(win, 1.5 * dt) + dt / 2).astype(bool)
            out["bkg_dc_excursion"] = hit.astype(float)

        # starguider: per telescope (m1/m2 only; LST has none) ; weather: site-wide
        starg = self._interp(dc_key, "starg", ("starg_dzd_deg", "starg_daz_deg"), t)
        out["starg_dev_deg"] = np.hypot(starg["starg_dzd_deg"], starg["starg_daz_deg"])
        wx_key = "m2" if "m2" in self.cond else next(iter(self.cond), None)
        if wx_key is not None:
            wx = self._interp(wx_key, "weather", ("wind_kmh", "wind_peak_kmh"), t)
            out["wind_kmh"], out["wind_peak_kmh"] = wx["wind_kmh"], wx["wind_peak_kmh"]
        return out

    def run_covariates(self, channel: str) -> Dict[str, float]:
        """Run-level wind / starguider covariates (90th percentiles)."""
        dc_key = self.night.telescope(channel).dc_key
        res = {}
        t_grid = np.linspace(self.t0, self.t1, 200)
        f = self.features(channel, t_grid)
        for k in ("starg_dev_deg", "wind_kmh", "wind_peak_kmh"):
            res[f"{k}_p90"] = float(np.nanpercentile(f[k], 90)) if np.isfinite(f[k]).any() else np.nan
        return res


def _near(t: np.ndarray, flagged: np.ndarray, pad: float) -> np.ndarray:
    """1.0 where ``t`` is within ``pad`` of any flagged time, else 0.0."""
    if not len(flagged):
        return np.zeros(len(t))
    flagged = np.sort(flagged)
    idx = np.searchsorted(flagged, t)
    prev = np.where(idx > 0, t - flagged[np.maximum(idx - 1, 0)], np.inf)
    nxt = np.where(idx < len(flagged), flagged[np.minimum(idx, len(flagged) - 1)] - t, np.inf)
    return (np.minimum(prev, nxt) <= pad).astype(float)


def _flags_from_features(f: Dict[str, np.ndarray], p90: Dict[str, float]) -> Dict[str, np.ndarray]:
    """Boolean flags used in the enrichment table."""
    def gt(x, thr):
        return np.where(np.isfinite(x) & np.isfinite(thr), x >= thr, False)
    return {
        "dc_excursion": f["dc_excursion"] > 0.5,
        "bkg_dc_excursion": f["bkg_dc_excursion"] > 0.5,
        "starg_dev_high": gt(f["starg_dev_deg"], p90.get("starg_dev_deg_p90", np.nan)),
        "wind_high": gt(f["wind_kmh"], p90.get("wind_kmh_p90", np.nan)),
        "wind_peak_high": gt(f["wind_peak_kmh"], p90.get("wind_peak_kmh_p90", np.nan)),
    }


# ─────────────────────────────── the pipeline ───────────────────────────────

@dataclass
class MimicResult:
    candidates: pd.DataFrame        # real + surrogate candidates (kind, k)
    events: pd.DataFrame            # coincidences (kind real / surrogate / timeslide)
    live: pd.DataFrame              # live time per (night, observation, correction, kind)
    tails: pd.DataFrame             # exceedance counts, real vs surrogate
    reference_flags: pd.DataFrame   # flag fractions at random reference times
    calibration: pd.DataFrame = field(default_factory=pd.DataFrame)   # SNR-trace std / MAD-sigma per run


def _correct(lc_full: LightCurve, name: str, night, dc_key: str, obs, cfg: MimicConfig):
    fn = _build_correction(name, {dc_key: night.dc_report} if night.dc_report else None, dc_key, 1.0)
    corrected = good_time_lightcurve(fn(lc_full), obs.good_time)
    return bin_lightcurve(corrected, cfg.bin_samples)


def run_mimic(log, nights=None, cfg: Optional[MimicConfig] = None, *, include_excluded: bool = False,
              progress: Optional[Callable[[str], None]] = None) -> MimicResult:
    """Run the full mimic measurement on the selected nights of the log."""
    cfg = cfg or MimicConfig()
    rng = np.random.default_rng(cfg.seed)
    say = progress or (lambda msg: None)
    cand_rows, event_rows, live_rows, tail_rows, ref_rows, calib_rows = [], [], [], [], [], []

    v = shadow_velocity_opposition(cfg.distance_au + EARTH_SEMIMAJOR_AXIS_AU)
    templates = search_templates(cfg, v)
    dur = _template_duration(templates)
    tol = coincidence_tolerance_s(ORM_ARRAY, ShadowKinematics(v), dur)

    for night in log.select(nights, include_excluded):
        for obs in night.stars.values():
            if not obs.has_cache:
                continue
            t0, t1 = obs.time_span()
            aux = AuxData(night, cfg, t0, t1)
            p90 = {ch: aux.run_covariates(ch) for ch in night.channels}
            lcs_full = {ch: obs.lightcurve(ch) for ch in night.channels}
            for corr in cfg.corrections:
                if corr in DC_CORRECTIONS and not night.dc_report:
                    continue
                say(f"{night.name} {obs.name} {corr}")
                lcs, real_rows = {}, {}
                try:
                    for ch in night.channels:
                        lcs[ch] = _correct(lcs_full[ch], corr, night, night.telescope(ch).dc_key, obs, cfg)
                except (ValueError, KeyError) as e:      # e.g. no DC data for this telescope
                    say(f"  skipped: {e}")
                    continue
                key = {"night": night.name, "observation": obs.name, "correction": corr}
                live = max(len(lc.time) * lc.dt for lc in lcs.values())
                live_rows.append({**key, "kind": "real", "live_s": live, "n_tel": len(lcs)})

                # ---- real data: candidates, aux features, tails
                real_tel = {}
                zs = {}
                for ch, lc in lcs.items():
                    rows, z, calib = search_lightcurve(lc, templates, cfg.snr_threshold, dur, cfg.sigma_window_s)
                    calib_rows.append({**key, "kind": "real", "channel": ch, **calib})
                    zs[ch] = z
                    real_tel[ch] = (rows, len(lc.time) * lc.dt)
                    feats = aux.features(ch, [r["time_s"] for r in rows]) if rows else None
                    for i, r in enumerate(rows):
                        flags = _flags_from_features({k: x[i:i + 1] for k, x in feats.items()}, p90[ch])
                        cand_rows.append({**key, "kind": "real", "k": -1, "channel": ch,
                                          "telescope": ORM_ARRAY.by_channel(ch).name, **r,
                                          **{k: float(x[i]) for k, x in feats.items()},
                                          **{f"flag_{k}": bool(x[0]) for k, x in flags.items()}})
                    # reference flag fractions at random times in the live span
                    tref = rng.uniform(lc.time[0], lc.time[-1], cfg.n_reference_times)
                    reff = aux.features(ch, tref)
                    for name, arr in _flags_from_features(reff, p90[ch]).items():
                        ref_rows.append({**key, "channel": ch, "flag": name,
                                         "n": int(np.isfinite(arr).sum()), "n_true": int(arr.sum())})

                res = _telescope_results(real_tel)
                if len(res) >= 2:
                    for e in _event_rows(res, tol):
                        event_rows.append({**key, "kind": "real", "k": -1, **e})
                    acc = accidental_coincidence_rate(res, tol, live, n_shifts=cfg.n_shifts,
                                                      min_shift_s=max(5.0, 3.0 * tol), seed=cfg.seed)
                    for snr in acc.event_min_snr:
                        event_rows.append({**key, "kind": "timeslide", "k": -1, "time_s": np.nan,
                                           "n_tel": len(res), "combined_snr": np.nan,
                                           "min_snr": float(snr), "max_chi2_reduced": np.nan,
                                           "channels": ""})
                    live_rows.append({**key, "kind": "timeslide", "live_s": live * cfg.n_shifts,
                                      "n_tel": len(res)})

                # ---- Gaussian surrogates, independent per channel
                surr_zs = {ch: [] for ch in lcs}
                live_rows.append({**key, "kind": "surrogate", "live_s": live * cfg.n_surrogates,
                                  "n_tel": len(lcs)})
                for k in range(cfg.n_surrogates):
                    surr_tel = {}
                    for ch, lc in lcs.items():
                        sur = phase_randomized_lightcurve(lc, rng)
                        rows, z, calib = search_lightcurve(sur, templates, cfg.snr_threshold, dur,
                                                           cfg.sigma_window_s)
                        calib_rows.append({**key, "kind": "surrogate", "channel": ch, **calib})
                        surr_zs[ch].append(z)
                        surr_tel[ch] = (rows, len(lc.time) * lc.dt)
                        for r in rows:
                            cand_rows.append({**key, "kind": "surrogate", "k": k, "channel": ch,
                                              "telescope": ORM_ARRAY.by_channel(ch).name, **r})
                    if len(surr_tel) >= 2:
                        for e in _event_rows(_telescope_results(surr_tel), tol):
                            event_rows.append({**key, "kind": "surrogate", "k": k, **e})

                # ---- tails (negative = dips, positive for reference), real vs surrogate
                for ch in lcs:
                    lc = lcs[ch]
                    n_eff = len(lc.time)
                    on = aux._dc_series(night.telescope(ch).dc_key, cfg.pixel_target)
                    for sign, label in ((-1, "dip"), (+1, "bump")):
                        for kσ in TAIL_SIGMAS:
                            idx = exceedances(zs[ch], kσ, sign)
                            n_real = len(idx)
                            n_sur = [len(exceedances(z, kσ, sign)) for z in surr_zs[ch]]
                            overlap = np.nan
                            if n_real and on is not None:
                                f = aux.features(ch, lc.time[idx])
                                overlap = float(np.mean((f["dc_excursion"] > 0.5) | (f["bkg_dc_excursion"] > 0.5)))
                            tail_rows.append({**key, "channel": ch, "tail": label, "sigma": kσ,
                                              "n_real": n_real, "n_surrogate_mean": float(np.mean(n_sur)),
                                              "n_surrogate_std": float(np.std(n_sur)),
                                              "n_white_gauss": float(n_eff * ndtr(-kσ)),
                                              "frac_overlap_dc_flag": overlap})

    return MimicResult(
        candidates=pd.DataFrame(cand_rows), events=pd.DataFrame(event_rows),
        live=pd.DataFrame(live_rows), tails=pd.DataFrame(tail_rows),
        reference_flags=pd.DataFrame(ref_rows), calibration=pd.DataFrame(calib_rows))


# ───────────────────────────────── tables ───────────────────────────────────

def far_table(result: MimicResult, snr_grid: Sequence[float] = (3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8)) -> pd.DataFrame:
    """
    False-alarm rate (events / hour) above each SNR, per correction, for
    ``scope`` = "single" (each telescope's candidates, per telescope-hour)
    and "coincident" (>= 2 telescopes; SNR = weakest member), for
    ``kind`` = real, surrogate and (coincident only) timeslide. rate_err is
    the Poisson error on the count.
    """
    rows = []
    cand, ev, live = result.candidates, result.events, result.live
    if cand.empty:
        return pd.DataFrame()
    for corr in sorted(cand["correction"].unique()):
        lv = live[live.correction == corr]
        for kind, scope in (("real", "single"), ("surrogate", "single"),
                            ("real", "coincident"), ("surrogate", "coincident"), ("timeslide", "coincident")):
            if scope == "single":
                sub = cand[(cand.correction == corr) & (cand.kind == kind)]
                # telescope-hours: live time of each run x number of telescopes
                lt = lv[lv.kind == kind]
                hours = float((lt.live_s * lt.n_tel).sum()) / 3600.0
                snr = sub["snr"].to_numpy()
            else:
                sub = ev[(ev.correction == corr) & (ev.kind == kind)] if not ev.empty else ev
                # only runs with >= 2 telescopes can have coincidences
                lt = lv[(lv.kind == kind) & (lv.n_tel >= 2)]
                hours = float(lt.live_s.sum()) / 3600.0
                snr = sub["min_snr"].to_numpy() if len(sub) else np.array([])
            if hours <= 0:
                continue
            for s in snr_grid:
                if scope == "coincident" and s < result_threshold(result):
                    continue
                n = int((snr >= s).sum())
                rows.append({"correction": corr, "scope": scope, "kind": kind, "snr_min": s,
                             "count": n, "hours": hours, "rate_per_hr": n / hours,
                             "rate_err": np.sqrt(max(n, 1)) / hours})
    return pd.DataFrame(rows)


def result_threshold(result: MimicResult) -> float:
    return float(result.candidates["snr"].min()) if not result.candidates.empty else 0.0


def tail_table(result: MimicResult) -> pd.DataFrame:
    """Exceedance counts summed over runs: real vs surrogate vs white Gaussian, with DC-flag overlap."""
    t = result.tails
    if t.empty:
        return t
    g = t.groupby(["correction", "tail", "sigma"])
    out = g.agg(n_real=("n_real", "sum"), n_surrogate=("n_surrogate_mean", "sum"),
                n_white_gauss=("n_white_gauss", "sum")).reset_index()
    out["excess_over_surrogate"] = out["n_real"] / out["n_surrogate"].replace(0, np.nan)
    w = t.assign(w=t["n_real"] * t["frac_overlap_dc_flag"]).groupby(["correction", "tail", "sigma"])["w"].sum()
    out["frac_real_with_dc_flag"] = (w.values / out["n_real"].replace(0, np.nan)).astype(float)
    return out


def flag_enrichment_table(result: MimicResult, snr_min: Sequence[float] = (3.0, 5.0)) -> pd.DataFrame:
    """
    For the real candidates above ``snr_min``: the fraction carrying each
    auxiliary flag, next to the fraction of random times carrying it
    (``base_fraction``); enrichment = fraction / base_fraction.
    """
    c, r = result.candidates, result.reference_flags
    c = c[c.kind == "real"] if not c.empty else c
    if c.empty or r.empty:
        return pd.DataFrame()
    base = (r.groupby(["correction", "flag"])[["n", "n_true"]].sum()
            .assign(base_fraction=lambda d: d.n_true / d.n.replace(0, np.nan))["base_fraction"])
    rows = []
    for corr in sorted(c.correction.unique()):
        for s in snr_min:
            sub = c[(c.correction == corr) & (c.snr >= s)]
            for col in [x for x in c.columns if x.startswith("flag_")]:
                name = col[5:]
                frac = float(sub[col].mean()) if len(sub) else np.nan
                b = float(base.get((corr, name), np.nan))
                rows.append({"correction": corr, "snr_min": s, "flag": name, "n_candidates": len(sub),
                             "fraction": frac, "base_fraction": b,
                             "enrichment": frac / b if b and np.isfinite(b) else np.nan})
    return pd.DataFrame(rows)
