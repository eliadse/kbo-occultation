# kbo_occultation/analysis/min_size.py
"""
Smallest detectable KBO of the array, predicted without injections.

The expected matched-filter SNR of a template in a measured noise spectrum is

    SNR = |G T| / sigma_mf(G T, |G|^2 S),

with T = I - 1 the template at the search cadence, S the noise PSD, G the (linear) noise
correction and sigma_mf = ``spectra.matched_filter_sigma`` (the scatter that the empirical
normalisation of ``matched_filter`` measures). Everything is in units of the star flux: the dip
is (1 - I) of the star signal and the noise is the fractional noise of the star signal, as in
``Instrument.noise_budget``.

Noise of a target star at any magnitude (``measure_noise_profiles``, ``noise_psd``):
the measured PSD of a real run, converted to star units (divided by the squared star fraction),
minus the white level of the noise model of that run, is the excess (correlated) noise;
it is assumed independent of magnitude and added to the white level that the noise model gives
for the target. ``theory_only`` drops the excess (ideal system).

Detection (``detection_probability``): each telescope's SNR is Gaussian around its expected
value (scatter ``snr_scatter``, 1 for the empirical normalisation) and must exceed a
threshold; the array needs ``min_ntel`` telescopes. Thresholds come from the real-data
false-alarm rate (``rate_curves``, ``coincidence_threshold``). The 50 % radius is ``r50``.

Caveats (to be calibrated by injections on the array, see min-size-plan.md):
- the excess noise inherits the model's NSB: if the measured NSB is above the model's,
  that difference sits in the excess;
- the rolling-median detrend of the search (20 template durations) is not modelled;
- real-data rate curves are the best of several templates per time, the prediction uses one.
"""

import json
from dataclasses import dataclass, field
from itertools import combinations
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfreqz
from scipy.special import ndtr

from .. import config
from ..config import (BandpassConfig, GridConfig, KBOConfig, NumericalConfig, ORM_ARRAY,
                      StarConfig)
from ..simulation.engine import OccultationEngine
from ..simulation.kinematics import ShadowKinematics, shadow_velocity_opposition
from ..simulation.physics import AU_m, mas_to_rad, teff_from_spectral_type
from .bands import SIGNAL_BAND_HZ
from .coincidence import coincidence_tolerance_s
from .spectra import lightcurve_psd, matched_filter_sigma

FRAME_S = config.standard_sampling * config.standard_sample_duration * 1e-9   # 524.288 us
SNR_FLOOR = 3.0                      # candidates are listed above this (noise_removal, mimic)
_TINY = 1e-300


# ───────────────────────────── targets ──────────────────────────────────────

@dataclass(frozen=True)
class TargetStar:
    name: str
    mag_B: float
    angular_radius_mas: float
    temperature_K: float = 10000.0
    spectral_type: Optional[str] = None


def star_radius_from_mag(mag_B: float, ref_mag: float = 7.0, ref_radius_mas: float = 0.03) -> float:
    """Angular radius [mas] at fixed surface brightness: r = r_ref * 10^(-0.2 (B - B_ref))."""
    return float(ref_radius_mas * 10.0 ** (-0.2 * (mag_B - ref_mag)))


def generic_targets(mags: Optional[Sequence[float]] = None, ref_mag: float = 7.0,
                    ref_radius_mas: float = 0.03, temperature_K: float = 10000.0) -> List[TargetStar]:
    """B = 7 ... 9 (step 0.25) stars whose radius shrinks with magnitude (star_radius_from_mag)."""
    mags = np.arange(7.0, 9.0 + 1e-9, 0.25) if mags is None else mags
    return [TargetStar(f"B{m:.2f}", float(m), star_radius_from_mag(m, ref_mag, ref_radius_mas), temperature_K)
            for m in mags]


CATALOGUE_COLUMNS = {"name": "name", "mag_B": "mag_B", "angular_radius_mas": "angular_radius_mas",
                     "temperature_K": "temperature_K", "spectral_type": "spectral_type"}


def load_target_catalogue(path: str, column_map: Optional[Mapping[str, str]] = None,
                          default_temperature_K: float = 10000.0, ref_mag: float = 7.0,
                          ref_radius_mas: float = 0.03) -> List[TargetStar]:
    """
    Targets from a CSV. ``column_map`` maps the fields (name, mag_B, angular_radius_mas,
    temperature_K, spectral_type) to the file's column names; absent columns are fine except
    mag_B. A missing radius comes from star_radius_from_mag, a missing temperature from the
    spectral type (physics.teff_from_spectral_type), else ``default_temperature_K``.
    """
    cols = {**CATALOGUE_COLUMNS, **(column_map or {})}
    df = pd.read_csv(path)
    if cols["mag_B"] not in df:
        raise ValueError(f"{path}: no magnitude column {cols['mag_B']!r} (set column_map['mag_B'])")

    def value(row, field_name):
        c = cols[field_name]
        return row[c] if c in df and pd.notna(row[c]) else None

    out = []
    for i, row in df.iterrows():
        mag = float(row[cols["mag_B"]])
        spt = value(row, "spectral_type")
        radius = value(row, "angular_radius_mas")
        temp = value(row, "temperature_K")
        if temp is None:
            temp = teff_from_spectral_type(spt) or default_temperature_K
        out.append(TargetStar(
            str(value(row, "name") or f"star{i}"), mag,
            float(radius) if radius is not None else star_radius_from_mag(mag, ref_mag, ref_radius_mas),
            float(temp), spt))
    return out


# ───────────────────────────── configuration ────────────────────────────────

@dataclass
class MinSizeConfig:
    distance_au: float = 39.0            # geocentric (Earth->KBO) at opposition
    heliocentric_au: float = 40.0        # sets the shadow velocity
    bin_samples: int = 10                # 5.24 ms, as in the search (mimic.MimicConfig)
    radii_m: np.ndarray = field(default_factory=lambda: np.geomspace(20.0, 1000.0, 30))
    correction: str = "whitened"         # raw | whitened | optimal | butter_<fc>
    coinc_far_per_hr: float = 1.0        # headline coincidence false-alarm rate
    extra_fars: Tuple[float, ...] = (0.1, 0.01)
    min_ntel: int = 2
    snr_scatter: float = 1.0             # sigma of the SNR of a telescope around its expectation
    snr_scale: float = 1.0               # factor on every expected SNR (to calibrate against injections)
    impact_parameter_m: float = 0.0
    excess_scaling: str = "star"         # excess noise fixed in star units ("star") or in total-signal units ("total")
    airmass: float = 1.2                 # used where a run has no airmass
    filter_type: Optional[str] = None
    x_half_m: float = 10000.0            # template half span (shared, as mimic.search_templates)
    dx_m: float = 5.0
    n_star_side: int = 64
    n_r_grid: int = 4000

    @property
    def bin_dt_s(self) -> float:
        return self.bin_samples * FRAME_S

    def velocity_mps(self) -> float:
        return shadow_velocity_opposition(self.heliocentric_au)

    def tolerance_s(self) -> float:
        """Coincidence window, as the real search (coincidence_tolerance_s, isotropic baseline)."""
        v = self.velocity_mps()
        return coincidence_tolerance_s(ORM_ARRAY, ShadowKinematics(v_rel_mps=v), 2.0 * self.x_half_m / v)


# ───────────────────────────── templates ────────────────────────────────────

def template_family(star: TargetStar, cfg: MinSizeConfig, instrument) -> Dict[float, np.ndarray]:
    """
    {radius_m: I - 1} at the bin cadence, for ``star`` seen through ``instrument.response()``
    (one engine, the Fresnel profiles are shared by the radii that repeat across calls).
    """
    v = cfg.velocity_mps()
    half = cfg.x_half_m
    grid = GridConfig(half, int(round(2.0 * half / cfg.dx_m)))
    r_star = star.angular_radius_mas * mas_to_rad * cfg.distance_au * AU_m
    r_max = 1.2 * np.hypot(half + r_star, r_star + abs(cfg.impact_parameter_m))
    band = BandpassConfig(*instrument.response_band())
    engine = OccultationEngine(StarConfig(star.temperature_K, star.angular_radius_mas), band, grid,
                               NumericalConfig(n_r_grid=cfg.n_r_grid, n_star_side=cfg.n_star_side),
                               response=instrument.response, r_max_m=r_max)
    t_fine = np.arange(-half / v, half / v, FRAME_S)
    t_fine = t_fine[:len(t_fine) // cfg.bin_samples * cfg.bin_samples]
    out = {}
    for radius in cfg.radii_m:
        x_m, intensity = engine.compute(KBOConfig(float(radius), cfg.distance_au, cfg.impact_parameter_m),
                                        star.angular_radius_mas)
        fine = np.interp(t_fine, x_m / v, intensity)
        out[float(radius)] = fine.reshape(-1, cfg.bin_samples).mean(axis=1) - 1.0
    return out


# ───────────────────────────── noise profiles ───────────────────────────────

@dataclass
class NoiseProfile:
    """
    Measured noise of one telescope on one star run, in star units (fractional noise of the
    star signal), plus what is needed to rebuild the noise model for another star under the
    same conditions.
    """
    night: str
    observation: str
    channel: str
    telescope: str
    mag_B: float
    teff_K: Optional[float]
    dt: float                      # cadence of the binned series [s]
    freq: np.ndarray               # Hz
    S_excess: np.ndarray           # 1/Hz, star units: measured PSD minus the model's white level, >= 0
    star_fraction: float           # S / (S + NSB) of the run (measured, or model if the NSB is missing)
    star_fraction_model: float
    sf_source: str                 # "measured" | "model"
    airmass: float
    psf_containment: float
    throughput: float
    k_B: Optional[float]
    clipped_fraction: float        # of the signal band where the measured PSD is below the model's white level
    live_s: float

    @property
    def fs(self) -> float:
        return 1.0 / self.dt

    _META = ("night", "observation", "channel", "telescope", "mag_B", "teff_K", "dt", "star_fraction",
             "star_fraction_model", "sf_source", "airmass", "psf_containment", "throughput", "k_B",
             "clipped_fraction", "live_s")


def model_budget(profile: NoiseProfile, mag_B: float, temperature_K: Optional[float],
                 focus_gain: float = 1.0) -> dict:
    """
    ``Instrument.noise_budget`` of a star of magnitude ``mag_B`` under the conditions of ``profile``
    (airmass with extinction, faulty-panel throughput, logged PSF, k_B of the site) at its cadence.
    ``focus_gain`` multiplies the star light of LST-1 (infinity focus); MAGIC is unchanged.
    """
    import astropy.units as u
    from .snr import instrument_for

    instr = instrument_for(profile.telescope, profile.k_B)
    g = focus_gain if instr.telescope_type == "LST" else 1.0
    return instr.noise_budget(mag_B, profile.dt * u.s, airmass=profile.airmass, extinction=True,
                              throughput_factor=profile.throughput, temperature_K=temperature_K,
                              psf_containment=profile.psf_containment * g)


def _star_fraction(budget: dict) -> float:
    return float(budget["star_photons"] / (budget["star_photons"] + budget["nsb_photons"]))


def measure_noise_profiles(log, nights=None, cfg: Optional[MinSizeConfig] = None, *,
                           include_excluded: bool = False, lcs=None, nperseg: int = 4096,
                           progress: Optional[Callable[[str], None]] = None) -> List[NoiseProfile]:
    """
    One NoiseProfile per (star run, channel) of the selected nights: the run is cropped to its
    good-time window, binned (``cfg.bin_samples``), and its Welch PSD converted to star units.
    ``nperseg`` is the segment length in binned samples (4096 = 21 s, 0.05 Hz resolution).
    The expensive step is loading the light curves (``lcs`` accepts the dict of
    ``snr.load_lightcurves``); use save_profiles to keep the result.
    """
    from ..observations.cc_reports import conditions_table
    from ..observations.photometry import good_time_lightcurve
    from .mimic import bin_lightcurve
    from .snr import instrument_for, load_lightcurves, logged_psf_containment, nsb_table
    from .snr_measurement import throughput_from_faulty_panels

    cfg = cfg or MinSizeConfig()
    say = progress or (lambda m: None)
    lcs = lcs if lcs is not None else load_lightcurves(log, nights, include_excluded=include_excluded)
    nsb = nsb_table(log, nights, lcs=lcs, units="fast", include_excluded=include_excluded)
    try:
        conditions = conditions_table(log, nights, include_excluded)
    except Exception as e:                                   # missing CC-report caches
        say(f"no conditions ({e}); airmass = {cfg.airmass}")
        conditions = pd.DataFrame()
    k_B = (getattr(log, "site", None) or {}).get("extinction_k", {}).get("B")
    lo, hi = SIGNAL_BAND_HZ

    profiles = []
    for night in log.select(nights, include_excluded):
        for obs in night.stars.values():
            if not obs.has_cache or obs.mag_B is None:
                continue
            airmass = conditions["airmass_mean"].get((night.name, obs.name), np.nan) \
                if "airmass_mean" in conditions else np.nan
            teff = teff_from_spectral_type(obs.spectral_type)
            for ch in night.channels:
                tel = night.telescope(ch)
                try:
                    lc = bin_lightcurve(good_time_lightcurve(lcs[(night.name, obs.name, ch)], obs.good_time),
                                        cfg.bin_samples)
                    freq, psd = lightcurve_psd(lc, nperseg)
                except (KeyError, ValueError) as e:
                    say(f"{night.name} {obs.name} {ch}: skipped ({e})")
                    continue
                profile = NoiseProfile(
                    night=night.name, observation=obs.name, channel=ch, telescope=tel.name, mag_B=obs.mag_B,
                    teff_K=teff, dt=float(lc.dt), freq=freq, S_excess=np.zeros_like(freq), star_fraction=np.nan,
                    star_fraction_model=np.nan, sf_source="model",
                    airmass=float(airmass) if np.isfinite(airmass) else cfg.airmass,
                    psf_containment=logged_psf_containment(instrument_for(tel.name, k_B),
                                                           (night.elog_notes or {}).get(tel.dc_key)),
                    throughput=throughput_from_faulty_panels(obs.faulty_mirrors.get(tel.dc_key), tel.dc_key),
                    k_B=k_B, clipped_fraction=np.nan, live_s=float(len(lc.time) * lc.dt))
                budget = model_budget(profile, obs.mag_B, teff)
                profile.star_fraction_model = _star_fraction(budget)
                sf = profile.star_fraction_model
                try:
                    nsb_fast = float(nsb.loc[(night.name, obs.name, tel.name), "nsb_default_fast"])
                    if np.isfinite(nsb_fast):
                        sf, profile.sf_source = 1.0 - nsb_fast / float(np.mean(lc.signal)), "measured"
                except KeyError:
                    pass
                if not 0.0 < sf <= 1.0:
                    say(f"{night.name} {obs.name} {ch}: star fraction {sf:.3g} outside (0, 1], using the model's")
                    sf, profile.sf_source = profile.star_fraction_model, "model"
                profile.star_fraction = float(sf)
                white = 2.0 * float(budget["sigma_total"]) ** 2 / profile.fs
                s_star = psd / sf ** 2
                profile.S_excess = np.maximum(s_star - white, 0.0)
                band = (freq >= lo) & (freq <= min(hi, freq[-1]))
                profile.clipped_fraction = float(np.mean(s_star[band] < white)) if band.any() else np.nan
                profiles.append(profile)
                say(f"{night.name} {obs.name} {ch}: sf={sf:.3f} ({profile.sf_source}), "
                    f"clipped {profile.clipped_fraction:.2f}")
    return profiles


def save_profiles(profiles: Sequence[NoiseProfile], path: str) -> None:
    meta = [{k: getattr(p, k) for k in NoiseProfile._META} for p in profiles]
    arrays = {}
    for i, p in enumerate(profiles):
        arrays[f"freq_{i}"], arrays[f"exc_{i}"] = p.freq, p.S_excess
    np.savez_compressed(path, meta=json.dumps(meta), **arrays)


def load_profiles(path: str) -> List[NoiseProfile]:
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z["meta"]))
        return [NoiseProfile(freq=z[f"freq_{i}"], S_excess=z[f"exc_{i}"], **m) for i, m in enumerate(meta)]


# ───────────────────────────── expected SNR ─────────────────────────────────

def correction_gain(name: str, freq: np.ndarray, psd: np.ndarray, fs: float) -> np.ndarray:
    """
    Amplitude gain G(f) of a linear noise correction on the rfft grid ``freq``:

    - ``raw``: 1;
    - ``whitened``: sqrt(min(1, S_white / S)), S_white the median PSD above min(100 Hz, 0.66 Nyquist)
      (noise_removal.whitening_response), 0 at f = 0;
    - ``optimal``: S^(-1/2), the best linear filter (upper bound);
    - ``butter_<fc>`` (e.g. ``butter_1``): the zero-phase 4th-order Butterworth high-pass of
      filtering.highpass_butterworth, |H|^2.
    """
    freq, psd = np.asarray(freq, float), np.asarray(psd, float)
    if name == "raw":
        return np.ones_like(freq)
    if name == "optimal":
        g = 1.0 / np.sqrt(np.maximum(psd, _TINY))
        return g / g.max()
    if name == "whitened":
        sel = freq >= min(100.0, 0.66 * fs / 2.0)
        s_white = float(np.median(psd[sel])) if sel.any() else float(np.min(psd[psd > 0]))
        h = np.sqrt(np.minimum(1.0, s_white / np.maximum(psd, _TINY)))
        h[0] = 0.0
        return h
    if name.startswith("butter_"):
        fc = float(name[len("butter_"):].removesuffix("Hz"))
        sos = butter(4, fc, btype="highpass", fs=fs, output="sos")
        _, h = sosfreqz(sos, worN=freq, fs=fs)
        return np.abs(h) ** 2
    raise ValueError(f"unknown correction {name!r}")


def expected_snr(template_residual: np.ndarray, freq: np.ndarray, psd: np.ndarray, fs: float,
                 gain: Optional[np.ndarray] = None) -> float:
    """
    Expected matched-filter SNR of the template (I - 1, sampled at 1 / fs) in the noise PSD ``psd``
    (one-sided, 1/Hz, on the rfft grid ``freq`` of n_fft = 2 (len(freq) - 1) points, as
    spectra.segmented_welch) after the correction ``gain`` (None = 1). For white noise of sigma s
    and no correction it is |T| / s.
    """
    n_fft = 2 * (len(freq) - 1)
    t = np.asarray(template_residual, dtype=float)
    if len(t) > n_fft:
        raise ValueError(f"template ({len(t)} samples) longer than the PSD's FFT ({n_fft})")
    g = np.ones_like(freq) if gain is None else np.asarray(gain, float)
    t_corr = np.fft.irfft(g * np.fft.rfft(t, n=n_fft), n=n_fft)
    sigma = matched_filter_sigma(freq, g ** 2 * np.asarray(psd, float), t_corr, fs)
    return float(np.linalg.norm(t_corr) / sigma)


def noise_psd(profile: NoiseProfile, star: TargetStar, cfg: MinSizeConfig, *, theory_only: bool = False,
              focus_gain: float = 1.0) -> np.ndarray:
    """
    Noise PSD (star units, 1/Hz, on profile.freq) of ``star`` in the conditions of ``profile``: the
    white level of the noise model plus the run's excess (not in ``theory_only``).
    """
    budget = model_budget(profile, star.mag_B, star.temperature_K, focus_gain)
    psd = np.full_like(profile.freq, 2.0 * float(budget["sigma_total"]) ** 2 / profile.fs)
    if theory_only:
        return psd
    excess = profile.S_excess
    if cfg.excess_scaling == "total":
        excess = excess * (profile.star_fraction_model / _star_fraction(budget)) ** 2
    elif cfg.excess_scaling != "star":
        raise ValueError(f"unknown excess_scaling {cfg.excess_scaling!r}")
    return psd + excess


def per_telescope_snr(star: TargetStar, profile: NoiseProfile, templates: Mapping[float, np.ndarray],
                      cfg: MinSizeConfig, *, theory_only: bool = False, focus_gain: float = 1.0) -> np.ndarray:
    """
    Expected SNR at every radius of cfg.radii_m times cfg.snr_scale (``theory_only``: white noise, no
    correction). Stationary Gaussian noise is assumed; against injections on the logged runs this
    is pessimistic (the real search normalises with a robust scatter, below the rms of
    intermittent red noise), hence ``snr_scale``.
    """
    if not np.isclose(profile.dt, cfg.bin_dt_s, rtol=1e-3):
        raise ValueError(f"profile cadence {profile.dt} s differs from the template's {cfg.bin_dt_s} s")
    psd = noise_psd(profile, star, cfg, theory_only=theory_only, focus_gain=focus_gain)
    gain = correction_gain("raw" if theory_only else cfg.correction, profile.freq, psd, profile.fs)
    return cfg.snr_scale * np.array([expected_snr(templates[float(r)], profile.freq, psd, profile.fs, gain)
                                     for r in cfg.radii_m])


# ───────────────────────────── thresholds from the real FAR ─────────────────

@dataclass
class RateCurve:
    """
    Rate of single-telescope candidates above an SNR (per hour) from the real-data candidate list:
    the empirical count up to the ``n_tail``-th largest candidate, an exponential fitted to the top
    ``n_tail`` beyond it (``extrapolated``). Fewer than 3 candidates: one per ``hours`` at the
    largest SNR with unit scale; none: zero rate (always extrapolated).
    """
    snr: np.ndarray
    hours: float
    n_tail: int = 20

    def __post_init__(self):
        a = np.sort(np.asarray(self.snr, dtype=float))[::-1]
        self.snr = a
        k = min(self.n_tail, len(a))
        if len(a) == 0:
            self._s_tail, self._amp, self._scale = np.inf, 0.0, 1.0
        elif k < 3:
            self._s_tail, self._amp, self._scale = a[0], 1.0 / self.hours, 1.0
        else:
            self._s_tail = a[k - 1]
            self._amp = (k - 1) / self.hours
            self._scale = max(float(np.mean(a[:k - 1] - a[k - 1])), 1e-3)

    def rate(self, s: float) -> float:
        if len(self.snr) == 0:
            return 0.0
        if s <= self._s_tail:
            return float(np.sum(self.snr > s)) / self.hours
        return self._amp * float(np.exp(-(s - self._s_tail) / self._scale))

    def extrapolated(self, s: float) -> bool:
        return len(self.snr) == 0 or s > self._s_tail


def csv_method_name(correction: str) -> str:
    """Name of a correction in outputs/noise_removal/real_candidates.csv."""
    if correction.startswith("butter_"):
        return f"butter_{float(correction[len('butter_'):].removesuffix('Hz')):g}Hz"
    return "whitened" if correction == "optimal" else correction


def rate_curves(candidates_csv: str, log, method: str, nights=None, *, include_excluded: bool = False,
                n_tail: int = 20) -> Dict[str, RateCurve]:
    """
    {channel: RateCurve} from ``real_candidates.csv`` (noise_removal.run_comparison) for one method
    (csv_method_name). Live time per channel is the good-time window of every cached star run of
    the selection that is in the CSV's nights (a night newer than the CSV would add live time
    without candidates).
    """
    df = pd.read_csv(candidates_csv)
    csv_nights = set(df["night"])
    df = df[df["method"] == method]
    hours: Dict[str, float] = {}
    runs = set()
    for night in log.select(nights, include_excluded):
        if night.name not in csv_nights:
            continue
        for obs in night.stars.values():
            if not obs.has_cache:
                continue
            t0, t1 = obs.time_span()
            span = (obs.good_time[1] - obs.good_time[0]) if obs.good_time else t1 - t0
            runs.add((night.name, obs.name))
            for ch in night.channels:
                hours[ch] = hours.get(ch, 0.0) + span / 3600.0
    df = df[[(n, o) in runs for n, o in zip(df["night"], df["observation"])]]
    return {ch: RateCurve(df.loc[df["channel"] == ch, "snr"].to_numpy(), h, n_tail) for ch, h in hours.items()}


def _solve_threshold(rate: Callable[[float], float], far: float, lo: float, hi: float = 100.0) -> float:
    if rate(lo) <= far:
        return lo
    if rate(hi) > far:
        return float("inf")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if rate(mid) > far else (lo, mid)
    return hi


def single_threshold(curve: RateCurve, far_per_hr: float, floor: float = SNR_FLOOR) -> Tuple[float, bool]:
    """(SNR threshold at which one telescope's real-data rate is ``far_per_hr``, extrapolated?)."""
    s = _solve_threshold(curve.rate, far_per_hr, floor)
    return s, curve.extrapolated(s)


def coincidence_rate(curves: Mapping[str, RateCurve], s: float, tolerance_s: float, min_ntel: int = 2) -> float:
    """
    Accidental rate [1/h] of coincidences of ``min_ntel`` of the telescopes at a common threshold ``s``:
    sum over telescope subsets of r_ref * prod_other (2 tau r_j), tau the coincidence window.
    """
    rates = {ch: c.rate(s) for ch, c in curves.items()}
    two_tau_h = 2.0 * tolerance_s / 3600.0
    total = 0.0
    for subset in combinations(rates, min_ntel):
        total += float(np.prod([rates[ch] for ch in subset])) * two_tau_h ** (min_ntel - 1)
    return total


def coincidence_threshold(curves: Mapping[str, RateCurve], tolerance_s: float, far_per_hr: float,
                          min_ntel: int = 2, floor: float = SNR_FLOOR) -> Tuple[float, bool]:
    """
    Common per-telescope SNR threshold at which the accidental ``min_ntel``-fold coincidence rate
    of independent streams is ``far_per_hr``, and whether any telescope's rate there is extrapolated.
    """
    if min_ntel > len(curves):
        raise ValueError(f"min_ntel={min_ntel} but only {len(curves)} telescopes")
    s = _solve_threshold(lambda x: coincidence_rate(curves, x, tolerance_s, min_ntel), far_per_hr, floor)
    return s, any(c.extrapolated(s) for c in curves.values())


# ───────────────────────────── detection probability and R50 ────────────────

def detection_probability(snr_by_tel: Mapping[str, np.ndarray], threshold, min_ntel: int,
                          snr_scatter: float = 1.0) -> np.ndarray:
    """
    P(at least ``min_ntel`` telescopes above their threshold) per radius, each SNR Gaussian around its
    expected value. ``threshold``: a number, or {telescope: number}.
    """
    probs = []
    for tel, snr in snr_by_tel.items():
        thr = threshold[tel] if isinstance(threshold, Mapping) else threshold
        probs.append(ndtr((np.asarray(snr, float) - thr) / snr_scatter))
    n = len(probs)
    dist = np.zeros((n + 1, len(probs[0])))
    dist[0] = 1.0
    for p in probs:                                            # Poisson-binomial recursion
        new = dist * (1.0 - p)
        new[1:] += dist[:-1] * p
        dist = new
    return dist[min_ntel:].sum(axis=0)


def r50(radii_m: Sequence[float], prob: Sequence[float], level: float = 0.5) -> Tuple[float, str]:
    """
    Smallest radius from which the probability stays >= ``level`` (log-interpolated), and a status:
    ``ok``; ``below_grid`` (already reached at the first radius, the value is that radius: an upper
    bound); ``above_grid`` (never reached, value nan).
    """
    radii, p = np.asarray(radii_m, float), np.asarray(prob, float)
    if p[-1] < level:
        return float("nan"), "above_grid"
    below = np.flatnonzero(p < level)
    if below.size == 0:
        return float(radii[0]), "below_grid"
    i = below[-1]
    frac = (level - p[i]) / (p[i + 1] - p[i])
    return float(np.exp(np.log(radii[i]) + frac * (np.log(radii[i + 1]) - np.log(radii[i])))), "ok"


# ───────────────────────────── the table ────────────────────────────────────

def min_size_table(targets: Sequence[TargetStar], profiles: Sequence[NoiseProfile], cfg: MinSizeConfig,
                   curves: Mapping[str, RateCurve], *, noise_models: Sequence[str] = ("measured", "theory"),
                   focus_gains: Sequence[float] = (1.0,),
                   progress: Optional[Callable[[str], None]] = None) -> pd.DataFrame:
    """
    R50 for every target x run (night, observation) x noise model x LST focus gain x telescope
    configuration (each telescope alone, ``k``of``n`` of the run's telescopes) x false-alarm rate
    (cfg.coinc_far_per_hr, then cfg.extra_fars).

    ``curves``: {channel: RateCurve} (rate_curves). The thresholds are those of the real data for both
    noise models, so ``theory`` isolates the effect of the noise on the SNR. Single telescopes use the
    same rate per telescope-hour as the coincidence FAR; arrays use coincidence_threshold. Rows:
    target, mag_B, star_mas, night, observation, noise, lst_focus, config, far_per_hr, threshold,
    thr_extrapolated, r50_m, r50_status, snr_100m (mean over the telescopes, interpolated).
    """
    from ..simulation.instruments import Instrument

    say = progress or (lambda m: None)
    fars = [cfg.coinc_far_per_hr, *cfg.extra_fars]
    tol = cfg.tolerance_s()
    runs: Dict[Tuple[str, str], List[NoiseProfile]] = {}
    for p in profiles:
        runs.setdefault((p.night, p.observation), []).append(p)
    instruments = {kind: Instrument(kind, filter_type=cfg.filter_type, site="LaPalma") for kind in ("MAGIC", "LST")}

    # thresholds depend only on the telescope set and the FAR
    thr_cache: Dict[tuple, Tuple[float, bool]] = {}

    def threshold(channels: Tuple[str, ...], k: int, far: float):
        key = (channels, k, far)
        if key not in thr_cache:
            sub = {ch: curves[ch] for ch in channels}
            thr_cache[key] = single_threshold(sub[channels[0]], far) if k == 1 else \
                coincidence_threshold(sub, tol, far, k)
        return thr_cache[key]

    rows = []
    for star in targets:
        templates = {kind: template_family(star, cfg, instr) for kind, instr in instruments.items()}
        say(f"{star.name}: B={star.mag_B:.2f}, r*={star.angular_radius_mas:.4f} mas")
        for (night, obs), group in runs.items():
            group = [p for p in group if p.channel in curves]
            for noise in noise_models:
                for gain in focus_gains:
                    if gain != 1.0 and not any(p.telescope.startswith("LST") for p in group):
                        continue
                    snr = {}
                    for p in group:
                        kind = "LST" if p.telescope.startswith("LST") else "MAGIC"
                        snr[p.channel] = per_telescope_snr(star, p, templates[kind], cfg,
                                                           theory_only=(noise == "theory"), focus_gain=gain)
                    chans = tuple(snr)
                    configs = [(p.telescope, (p.channel,), 1) for p in group]
                    configs += [(f"{k}of{len(chans)}", chans, k) for k in range(2, len(chans) + 1)]
                    for label, sel, k in configs:
                        mean_snr = np.mean([snr[c] for c in sel], axis=0)
                        for far in fars:
                            thr, extrap = threshold(sel, k, far)
                            prob = detection_probability({c: snr[c] for c in sel}, thr, k, cfg.snr_scatter)
                            value, status = r50(cfg.radii_m, prob)
                            rows.append({
                                "target": star.name, "mag_B": star.mag_B, "star_mas": star.angular_radius_mas,
                                "night": night, "observation": obs, "noise": noise, "lst_focus": gain,
                                "config": label, "far_per_hr": far, "threshold": thr,
                                "thr_extrapolated": extrap, "r50_m": value, "r50_status": status,
                                "snr_100m": float(np.interp(np.log(100.0), np.log(cfg.radii_m), mean_snr))})
    return pd.DataFrame(rows)


def summarize(table: pd.DataFrame) -> pd.DataFrame:
    """Median / min / max R50 over the runs, per target, noise, focus, configuration and FAR."""
    keys = ["target", "mag_B", "star_mas", "noise", "lst_focus", "config", "far_per_hr"]
    g = table.groupby(keys)
    return g.agg(r50_median=("r50_m", "median"), r50_min=("r50_m", "min"), r50_max=("r50_m", "max"),
                 n_runs=("r50_m", "size"), n_unreached=("r50_status", lambda s: int((s == "above_grid").sum())),
                 n_bound=("r50_status", lambda s: int((s == "below_grid").sum()))).reset_index()
