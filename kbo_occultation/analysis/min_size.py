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

Calibration (``calibrate_run`` ... ``fit_snr_calibration``, driven by
examples/studies/min_size_calibration.py): the same event is injected into all channels of a logged
run, searched as the real data are, and the measured efficiency / peak SNR compared with the
prediction above; the fitted (scale, scatter) per correction go to ``SNR_CALIBRATION``.

Caveats (to be calibrated by injections on the array, see min-size-plan.md):
- the excess noise inherits the model's NSB: if the measured NSB is above the model's,
  that difference sits in the excess;
- the rolling-median detrend of the search (20 template durations) is not modelled;
- real-data rate curves are the best of several templates per time, the prediction uses one.
"""

import json
import zlib
from dataclasses import dataclass, field, replace
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
from .spectra import _trapz, lightcurve_psd, matched_filter_sigma

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
    snr_scale: float = 1.0               # factor on the measured-noise SNR; ``calibrated()`` takes SNR_CALIBRATION's
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

    def calibrated(self) -> "MinSizeConfig":
        """Copy with ``snr_scale`` / ``snr_scatter`` of this correction from ``SNR_CALIBRATION`` (unchanged if absent)."""
        scale, scatter = SNR_CALIBRATION.get(self.correction, (self.snr_scale, self.snr_scatter))
        return replace(self, snr_scale=scale, snr_scatter=scatter)

    def velocity_mps(self) -> float:
        return shadow_velocity_opposition(self.heliocentric_au)

    def tolerance_s(self) -> float:
        """Coincidence window, as the real search (coincidence_tolerance_s, isotropic baseline)."""
        v = self.velocity_mps()
        return coincidence_tolerance_s(ORM_ARRAY, ShadowKinematics(v_rel_mps=v), 2.0 * self.x_half_m / v)


# ───────────────────────────── templates ────────────────────────────────────

def _engine(star: TargetStar, cfg: MinSizeConfig, instrument) -> OccultationEngine:
    """Occultation engine of ``star`` seen through ``instrument.response()`` on the shared template grid."""
    half = cfg.x_half_m
    grid = GridConfig(half, int(round(2.0 * half / cfg.dx_m)))
    r_star = star.angular_radius_mas * mas_to_rad * cfg.distance_au * AU_m
    r_max = 1.2 * np.hypot(half + r_star, r_star + abs(cfg.impact_parameter_m))
    band = BandpassConfig(*instrument.response_band())
    return OccultationEngine(StarConfig(star.temperature_K, star.angular_radius_mas), band, grid,
                             NumericalConfig(n_r_grid=cfg.n_r_grid, n_star_side=cfg.n_star_side),
                             response=instrument.response, r_max_m=r_max)


def injection_patterns(star: TargetStar, cfg: MinSizeConfig, instrument,
                       radii: Optional[Sequence[float]] = None) -> Dict[float, Tuple[np.ndarray, np.ndarray]]:
    """{radius_m: (x_m, I)}: the Fresnel pattern of ``star`` through ``instrument`` at full resolution, for injection."""
    engine = _engine(star, cfg, instrument)
    out = {}
    for radius in (cfg.radii_m if radii is None else radii):
        out[float(radius)] = engine.compute(KBOConfig(float(radius), cfg.distance_au, cfg.impact_parameter_m),
                                            star.angular_radius_mas)
    return out


def template_family(star: TargetStar, cfg: MinSizeConfig, instrument) -> Dict[float, np.ndarray]:
    """
    {radius_m: I - 1} at the bin cadence, for ``star`` seen through ``instrument.response()``
    (one engine, the Fresnel profiles are shared by the radii that repeat across calls).
    """
    v = cfg.velocity_mps()
    half = cfg.x_half_m
    t_fine = np.arange(-half / v, half / v, FRAME_S)
    t_fine = t_fine[:len(t_fine) // cfg.bin_samples * cfg.bin_samples]
    out = {}
    for radius, (x_m, intensity) in injection_patterns(star, cfg, instrument).items():
        fine = np.interp(t_fine, x_m / v, intensity)
        out[radius] = fine.reshape(-1, cfg.bin_samples).mean(axis=1) - 1.0
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
    Expected SNR at every radius of cfg.radii_m. With measured noise it is multiplied by
    cfg.snr_scale (see ``SNR_CALIBRATION``: the real search normalises with a robust scatter, which on
    intermittent red noise is below the rms that the PSD gives); ``theory_only`` (white noise, no
    correction) is not scaled, the scale was fitted on real data. Stationary Gaussian noise is assumed.
    """
    if not np.isclose(profile.dt, cfg.bin_dt_s, rtol=1e-3):
        raise ValueError(f"profile cadence {profile.dt} s differs from the template's {cfg.bin_dt_s} s")
    psd = noise_psd(profile, star, cfg, theory_only=theory_only, focus_gain=focus_gain)
    gain = correction_gain("raw" if theory_only else cfg.correction, profile.freq, psd, profile.fs)
    scale = 1.0 if theory_only else cfg.snr_scale
    return scale * np.array([expected_snr(templates[float(r)], profile.freq, psd, profile.fs, gain)
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


# ───────────────────────────── calibration against injections ───────────────

# correction -> (snr_scale, snr_scatter) fitted on array injections (examples/studies/min_size_calibration.py,
# outputs/min_size/calibration/snr_calibration.csv). 2026-11-12: 2026-03-08 + 2026-08-19, 7 runs x 3 channels,
# 50 trials x 12 radii (50-400 m), measured star fraction, FAR 1/h. Median over runs of (measured peak SNR) /
# (step-1 prediction), per-run range 0.67-1.65 (whitened), 0.89-1.58 (butter_1), 0.87-6.5 (raw); scatter is the
# absolute 1.4826 MAD of the peak SNR below SNR 20. The likelihood fit on the detections gives 0.96 / 1.15 / 1.37.
SNR_CALIBRATION: Dict[str, Tuple[float, float]] = {"whitened": (0.912, 1.234), "raw": (1.539, 2.250),
                                                   "butter_1": (1.082, 1.463)}

CALIBRATION_RADII_M = (50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 400.0)
_RUN = ["night", "observation", "correction", "sf_mode"]


def observed_target(obs, star_mas: float = 0.03) -> TargetStar:
    """The star of a logged run as a TargetStar (the injections and the search bank use 0.03 mas)."""
    teff = teff_from_spectral_type(obs.spectral_type) if obs.spectral_type else None
    return TargetStar(obs.name, float(obs.mag_B), star_mas, float(teff or 10000.0), obs.spectral_type)


def calibration_runs(log, nights, profiles: Sequence[NoiseProfile], *, include_excluded: bool = False,
                     min_channels: int = 2):
    """[(night, observation, [NoiseProfile per channel])] of the cached runs that have >= ``min_channels`` profiles."""
    by_run: Dict[Tuple[str, str], List[NoiseProfile]] = {}
    for p in profiles:
        by_run.setdefault((p.night, p.observation), []).append(p)
    out = []
    for night in log.select(nights, include_excluded):
        for obs in night.stars.values():
            group = by_run.get((night.name, obs.name), [])
            if obs.has_cache and len(group) >= min_channels:
                out.append((night, obs, sorted(group, key=lambda p: p.channel)))
    return out


def _is_lst(p: NoiseProfile) -> bool:
    return p.telescope.startswith("LST")


def calibrate_run(night, obs, group: Sequence[NoiseProfile], cfg: MinSizeConfig, ccfg, correction: str,
                  radii: Sequence[float] = CALIBRATION_RADII_M, *, sf_mode: str = "measured",
                  n_trials: Optional[int] = None, seed: int = 0, bank=None,
                  patterns: Optional[dict] = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Array injections of one run and one correction (``noise_removal.injection_efficiency_array``).

    The search is the one that wrote real_candidates.csv: the method's reshaped bank of
    ``mimic.search_templates(ccfg.mimic, ...)``. The injected pattern is the observed star through the
    telescope's own response (``injection_patterns``), diluted by the star fraction of each channel:
    ``sf_mode="measured"`` is the profile's (what the prediction uses), ``"model"`` the noise
    model's. ``patterns`` is a {kind: {R: (x, I)}} cache that is filled in (pass the same dict for the
    corrections of one run and the same ``radii``). Returns (trials, candidates, run_calib);
    run_calib has the full-run search statistics per channel (SNR-trace std and MAD-sigma, number of
    candidates) to compare with real_candidates.csv.
    """
    from ..observations.photometry import good_time_lightcurve
    from ..simulation.instruments import Instrument
    from ..simulation.kinematics import EARTH_SEMIMAJOR_AXIS_AU
    from .mimic import _template_duration, bin_lightcurve, search_lightcurve, search_templates
    from .noise_removal import (build_methods, injection_efficiency_array, reshaped_templates,
                                shadow_direction_enu)

    if ccfg.mimic.bin_samples != cfg.bin_samples:
        raise ValueError("ComparisonConfig.mimic.bin_samples and MinSizeConfig.bin_samples differ")
    if sf_mode not in ("measured", "model"):
        raise ValueError(f"unknown sf_mode {sf_mode!r}")
    if n_trials is not None:
        ccfg = replace(ccfg, n_trials=n_trials)
    name = csv_method_name(correction)
    v_bank = shadow_velocity_opposition(ccfg.mimic.distance_au + EARTH_SEMIMAJOR_AXIS_AU)
    templates_raw = bank if bank is not None else search_templates(ccfg.mimic, v_bank)
    dur = _template_duration(templates_raw)
    channels = [p.channel for p in group]
    lcs = {ch: obs.lightcurve(ch) for ch in channels}
    key = {"night": night.name, "observation": obs.name, "correction": correction, "sf_mode": sf_mode}

    methods, bank_m, corrected, calib_rows = {}, {}, {}, []
    for p in group:
        ch = p.channel
        methods[ch] = build_methods(lcs[ch], night.dc_report, night.telescope(ch).dc_key, [name])[name]
        bank_m[ch] = reshaped_templates(lcs[ch], methods[ch], templates_raw, v_bank, ccfg)
        corrected[ch] = methods[ch].correct(lcs[ch])
        real = bin_lightcurve(good_time_lightcurve(corrected[ch], obs.good_time), ccfg.mimic.bin_samples)
        rows, _, c = search_lightcurve(real, bank_m[ch], ccfg.snr_floor, dur, ccfg.mimic.sigma_window_s)
        calib_rows.append({**key, "channel": ch, "telescope": p.telescope, "n_candidates": len(rows),
                           "n_ge5": sum(r["snr"] >= 5.0 for r in rows), "snr_std": c["snr_std"],
                           "snr_mad_sigma": c["snr_mad_sigma"], "live_s": len(real.time) * real.dt})

    star = observed_target(obs)
    patterns = {} if patterns is None else patterns
    for kind in {"LST" if _is_lst(p) else "MAGIC" for p in group}:
        if kind not in patterns:
            patterns[kind] = injection_patterns(star, cfg, Instrument(kind, filter_type=cfg.filter_type,
                                                                      site="LaPalma"), radii)
    ch_patterns = {p.channel: patterns["LST" if _is_lst(p) else "MAGIC"] for p in group}
    fractions = {p.channel: (p.star_fraction if sf_mode == "measured" else p.star_fraction_model) for p in group}

    t_start = obs.time_span()[0]
    t_mid = t_start + (float(np.mean(obs.good_time)) if obs.good_time else 0.0)
    try:
        direction = shadow_direction_enu(obs.ra_deg, obs.dec_deg, t_mid, cfg.distance_au)
    except Exception:                                         # no ephemeris / astropy data: isotropic
        direction = None
    rng = np.random.default_rng([seed, zlib.crc32(f"{night.name}/{obs.name}/{name}/{sf_mode}".encode())])
    trial_rows, cand_rows = injection_efficiency_array(
        lcs, methods, ch_patterns, bank_m, dur, cfg.velocity_mps(), fractions, obs.good_time, ccfg, rng,
        list(radii), direction_enu=direction, corrected_full=corrected)
    return (pd.DataFrame([{**key, **r} for r in trial_rows]), pd.DataFrame([{**key, **r} for r in cand_rows]),
            pd.DataFrame(calib_rows))


# ---- false-alarm check of the thresholds

def timeslide_far(candidates_csv: str, log, method: str, nights=None, tolerance_s: Optional[float] = None,
                  snr_grid: Sequence[float] = tuple(np.arange(3.0, 12.01, 0.25)), n_shifts: int = 2000,
                  min_ntel: int = 2, seed: int = 0, curves: Optional[Mapping[str, RateCurve]] = None,
                  include_excluded: bool = False) -> pd.DataFrame:
    """
    Coincidence false-alarm rate of the real candidate lists (``real_candidates.csv`` rows of ``method``)
    on the runs with >= ``min_ntel`` channels, measured with time slides
    (``coincidence.accidental_coincidence_rate``), against the zero-lag count and the independent-stream
    estimate ``coincidence_rate(curves, s, tolerance_s)``. The SNR of an event is its weakest member's, so
    three-fold chance events (third order) are counted at their lowest SNR.

    Columns: snr, slide_count, slide_hours, slide_rate_per_hr, zero_lag_count, live_hours,
    zero_lag_rate_per_hr, analytic_rate_per_hr (nan without ``curves``).
    """
    from .coincidence import TelescopeSearchResult, accidental_coincidence_rate, match_coincidences
    from .matched_filter import Candidate

    tol = tolerance_s if tolerance_s is not None else MinSizeConfig().tolerance_s()
    df = pd.read_csv(candidates_csv)
    df = df[df["method"] == method]
    slide_snr, zero_snr = [], []
    slide_h = live_h = 0.0
    for night in log.select(nights, include_excluded):
        if len(night.channels) < min_ntel:
            continue
        for obs in night.stars.values():
            if not obs.has_cache:
                continue
            t0, t1 = obs.time_span()
            a, b = (t0 + obs.good_time[0], t0 + obs.good_time[1]) if obs.good_time else (t0, t1)
            live, min_shift = b - a, max(5.0, 3.0 * tol)
            if live <= 2.0 * min_shift:
                continue
            sub = df[(df["night"] == night.name) & (df["observation"] == obs.name)]
            results = [TelescopeSearchResult(
                ch, night.telescope(ch).name,
                [Candidate(time_s=r.time_s, snr=r.snr, chi2_reduced=r.chi2_reduced)
                 for r in sub[sub["channel"] == ch].itertuples()], 0.0, live) for ch in night.channels]
            acc = accidental_coincidence_rate(results, tol, live, n_shifts=n_shifts, min_shift_s=min_shift,
                                              min_ntel=min_ntel, seed=seed, t_start_s=a)
            slide_snr.extend(acc.event_min_snr)
            slide_h += n_shifts * live / 3600.0
            zero_snr.extend(min(c.snr for c in e.members.values())
                            for e in match_coincidences(results, tol, min_ntel=min_ntel))
            live_h += live / 3600.0
    slide_snr, zero_snr = np.asarray(slide_snr, float), np.asarray(zero_snr, float)
    rows = []
    for snr in snr_grid:
        ns, nz = int((slide_snr >= snr).sum()), int((zero_snr >= snr).sum())
        rows.append({"snr": float(snr), "slide_count": ns, "slide_hours": slide_h,
                     "slide_rate_per_hr": ns / slide_h if slide_h else np.nan, "zero_lag_count": nz,
                     "live_hours": live_h, "zero_lag_rate_per_hr": nz / live_h if live_h else np.nan,
                     "analytic_rate_per_hr": coincidence_rate(curves, float(snr), tol, min_ntel) if curves else np.nan})
    return pd.DataFrame(rows)


def thresholds_table(curves: Mapping[str, RateCurve], tolerance_s: float, fars: Sequence[float],
                     far_check: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Per-telescope SNR thresholds for each false-alarm rate: every channel alone (``single_threshold``) and
    k of n (``coincidence_threshold``, common value), flagged when extrapolated beyond the measured tail.
    On the 2-of-n rows ``timeslide_threshold`` is the lowest SNR whose time-slide rate (``far_check``) is
    <= the FAR; ``timeslide_resolved`` says whether the slides could have seen one event at that rate.
    """
    n = len(curves)
    rows = []
    for far in fars:
        for ch, curve in curves.items():
            thr, ext = single_threshold(curve, far)
            rows.append({"far_per_hr": far, "config": ch, "k": 1, "threshold": thr, "thr_extrapolated": ext,
                         "timeslide_threshold": np.nan, "timeslide_resolved": np.nan})
        for k in range(2, n + 1):
            thr, ext = coincidence_threshold(curves, tolerance_s, far, k)
            ts, resolved = np.nan, np.nan
            if k == 2 and far_check is not None and len(far_check):
                ok = far_check[far_check["slide_rate_per_hr"] <= far]
                ts = float(ok["snr"].min()) if len(ok) else np.nan
                resolved = bool(far * float(far_check["slide_hours"].iloc[0]) >= 1.0)
            rows.append({"far_per_hr": far, "config": f"{k}of{n}", "k": k, "threshold": thr,
                         "thr_extrapolated": ext, "timeslide_threshold": ts, "timeslide_resolved": resolved})
    return pd.DataFrame(rows)


# ---- measured and predicted efficiency

def _binomial_err(hits: int, n: int) -> float:
    p = min(max(hits / n, 0.5 / n), 1.0 - 0.5 / n) if n else np.nan
    return float(np.sqrt(p * (1.0 - p) / n)) if n else np.nan


def measured_efficiency(trials: pd.DataFrame, cands: pd.DataFrame, thresholds: pd.DataFrame,
                        tolerance_s: float) -> pd.DataFrame:
    """
    Fraction of the injections detected, per run, correction, radius, configuration and false-alarm rate.

    Single channel: it has a candidate (within the search window of ``injection_efficiency_array``) above its
    threshold. k of n: ``coincidence.match_coincidences(min_ntel=k)`` of the candidates above the common
    threshold gives an event within ``tolerance_s`` of the injection time. Rows: the run keys, radius_m,
    config, far_per_hr, threshold, n, hits, eff, eff_err.
    """
    from .coincidence import TelescopeSearchResult, match_coincidences
    from .matched_filter import Candidate

    if trials.empty:
        return pd.DataFrame()
    thr = {(r.far_per_hr, r.config): r.threshold for r in thresholds.itertuples()}
    fars = sorted({f for f, _ in thr})
    grp_cols = _RUN + ["radius_m", "trial"]
    t0s = trials.groupby(grp_cols)["t0_s"].mean()
    chans = trials.groupby(_RUN)["channel"].agg(lambda c: sorted(set(c)))
    cand_by = {k: (g["channel"].to_numpy(), g["snr"].to_numpy(), g["time_s"].to_numpy(), g["chi2_reduced"].to_numpy())
               for k, g in cands.groupby(grp_cols)} if len(cands) else {}
    empty = (np.array([], dtype=object), np.array([]), np.array([]), np.array([]))

    hits: Dict[tuple, int] = {}
    count: Dict[tuple, int] = {}
    for key, t0 in t0s.items():
        run, radius = key[:4], key[4]
        channels = chans[run]
        ch_arr, snr_arr, t_arr, chi_arr = cand_by.get(key, empty)
        for far in fars:
            configs = [(c, 1) for c in channels] + [(f"{k}of{len(channels)}", k) for k in range(2, len(channels) + 1)]
            for config, k in configs:
                if (far, config) not in thr:
                    continue
                s = thr[(far, config)]
                if k == 1:
                    ok = bool(np.any((ch_arr == config) & (snr_arr >= s)))
                else:
                    res = []
                    for ch in channels:
                        m = (ch_arr == ch) & (snr_arr >= s)
                        res.append(TelescopeSearchResult(
                            ch, ch, [Candidate(time_s=float(t), snr=float(q), chi2_reduced=float(c))
                                     for t, q, c in zip(t_arr[m], snr_arr[m], chi_arr[m])], 0.0, 0.0))
                    ok = any(abs(e.time_s - t0) <= tolerance_s for e in match_coincidences(res, tolerance_s, k))
                kk = (*run, radius, config, far)
                count[kk] = count.get(kk, 0) + 1
                hits[kk] = hits.get(kk, 0) + int(ok)
    rows = [{**dict(zip(_RUN, kk[:4])), "radius_m": kk[4], "config": kk[5], "far_per_hr": kk[6],
             "threshold": thr[(kk[6], kk[5])], "n": n, "hits": hits[kk], "eff": hits[kk] / n,
             "eff_err": _binomial_err(hits[kk], n)} for kk, n in count.items()]
    return pd.DataFrame(rows)


def predicted_snr(runs, cfg: MinSizeConfig, corrections: Sequence[str],
                  radii: Sequence[float] = CALIBRATION_RADII_M) -> pd.DataFrame:
    """
    Step-1 expected SNR (``per_telescope_snr``, snr_scale 1) of the observed star of every run at the
    injection radii, per correction and channel: night, observation, correction, channel, telescope,
    radius_m, snr_pred. ``runs``: ``calibration_runs``.
    """
    from ..simulation.instruments import Instrument

    cfg_r = replace(cfg, radii_m=np.asarray(radii, float), snr_scale=1.0)
    instruments = {k: Instrument(k, filter_type=cfg.filter_type, site="LaPalma") for k in ("MAGIC", "LST")}
    rows = []
    for night, obs, group in runs:
        star = observed_target(obs)
        kinds = {"LST" if _is_lst(p) else "MAGIC" for p in group}
        templates = {k: template_family(star, cfg_r, instruments[k]) for k in kinds}
        for corr in corrections:
            cfg_c = replace(cfg_r, correction=corr)
            for p in group:
                snr = per_telescope_snr(star, p, templates["LST" if _is_lst(p) else "MAGIC"], cfg_c)
                rows.extend({"night": night.name, "observation": obs.name, "correction": corr,
                             "channel": p.channel, "telescope": p.telescope, "radius_m": float(r),
                             "snr_pred": float(v)} for r, v in zip(radii, snr))
    return pd.DataFrame(rows)


def predicted_efficiency(snr_pred: pd.DataFrame, thresholds: pd.DataFrame,
                         scales: Optional[Mapping[str, Tuple[float, float]]] = None) -> pd.DataFrame:
    """
    Detection probability from ``predicted_snr`` at the same thresholds as ``measured_efficiency``, with the
    expected SNRs multiplied by ``scales[correction][0]`` and a scatter ``scales[correction][1]`` (default
    1, 1). Rows: night, observation, correction, radius_m, config, far_per_hr, threshold, eff_pred.
    """
    scales = scales or {}
    thr = {(r.far_per_hr, r.config): r.threshold for r in thresholds.itertuples()}
    fars = sorted({f for f, _ in thr})
    rows = []
    for (night, obs, corr), g in snr_pred.groupby(["night", "observation", "correction"]):
        scale, scatter = scales.get(corr, (1.0, 1.0))
        radii = np.sort(g["radius_m"].unique())
        snr = {ch: gg.sort_values("radius_m")["snr_pred"].to_numpy() * scale for ch, gg in g.groupby("channel")}
        n = len(snr)
        for far in fars:
            for ch, v in snr.items():
                if (far, ch) in thr:
                    p = ndtr((v - thr[(far, ch)]) / scatter)
                    rows.extend({"night": night, "observation": obs, "correction": corr, "radius_m": float(r),
                                 "config": ch, "far_per_hr": far, "threshold": thr[(far, ch)], "eff_pred": float(x)}
                                for r, x in zip(radii, p))
            for k in range(2, n + 1):
                label = f"{k}of{n}"
                if (far, label) in thr:
                    p = detection_probability(snr, thr[(far, label)], k, scatter)
                    rows.extend({"night": night, "observation": obs, "correction": corr, "radius_m": float(r),
                                 "config": label, "far_per_hr": far, "threshold": thr[(far, label)],
                                 "eff_pred": float(x)} for r, x in zip(radii, p))
    return pd.DataFrame(rows)


# ---- calibration fit

def fit_snr_calibration(trials: pd.DataFrame, snr_pred: pd.DataFrame, run_calib: Optional[pd.DataFrame] = None,
                        snr_range: Tuple[float, float] = (5.0, 50.0), scatter_max: float = 20.0
                        ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Scale and scatter of the measured peak SNR against the prediction. Per (radius, run, channel) the cell's
    median peak SNR over the trials divided by ``snr_pred`` is a scale estimate; cells with a median in
    ``snr_range`` are used (below it the maximum over the search window is biased high by noise, above it the
    SNR is saturated by shape mismatch). The run's scale is the median over its radii; ``scatter`` is the
    median 1.4826 MAD of peak SNR within the cells with a median <= ``scatter_max`` (absolute SNR units, as
    ``MinSizeConfig.snr_scatter``).

    Returns (per_run, summary): per_run rows are (run keys, channel, telescope, scale, scatter, n_cells,
    snr_std, snr_mad_sigma) with the full-run SNR-trace statistics of ``run_calib`` (the robust-scatter
    hypothesis: scale ~ snr_std); summary is the median over runs per correction, sf_mode and telescope kind.
    """
    t = trials.dropna(subset=["peak_snr"])
    mad = lambda x: 1.4826 * float(np.median(np.abs(x - np.median(x))))                   # noqa: E731
    cells = (t.groupby(_RUN + ["channel", "radius_m"])["peak_snr"]
             .agg(median="median", mad=mad, n="count").reset_index()
             .merge(snr_pred, on=["night", "observation", "correction", "channel", "radius_m"]))
    cells["ratio"] = cells["median"] / cells["snr_pred"]
    use = cells[(cells["median"] >= snr_range[0]) & (cells["median"] <= snr_range[1])]
    keys = _RUN + ["channel", "telescope"]
    per_run = use.groupby(keys).agg(scale=("ratio", "median"), n_cells=("ratio", "size")).reset_index()
    low = use[use["median"] <= scatter_max].groupby(keys)["mad"].median().rename("scatter").reset_index()
    per_run = per_run.merge(low, on=keys, how="left")
    if run_calib is not None and len(run_calib):
        per_run = per_run.merge(run_calib[["night", "observation", "correction", "sf_mode", "channel",
                                           "snr_std", "snr_mad_sigma"]].drop_duplicates(),
                                on=_RUN + ["channel"], how="left")
    per_run["kind"] = np.where(per_run["telescope"].str.startswith("LST"), "LST", "MAGIC")
    summary = (per_run.groupby(["correction", "sf_mode", "kind"])
               .agg(scale=("scale", "median"), scale_min=("scale", "min"), scale_max=("scale", "max"),
                    scatter=("scatter", "median"), n_runs=("scale", "size"),
                    snr_std=("snr_std", "median") if "snr_std" in per_run else ("scale", "size"))
               .reset_index())
    pooled = (per_run.groupby(["correction", "sf_mode"])
              .agg(scale=("scale", "median"), scale_min=("scale", "min"), scale_max=("scale", "max"),
                   scatter=("scatter", "median"), n_runs=("scale", "size"),
                   snr_std=("snr_std", "median") if "snr_std" in per_run else ("scale", "size"))
              .reset_index().assign(kind="all"))
    return per_run, pd.concat([summary, pooled], ignore_index=True)


def scale_by_deviance(eff_meas: pd.DataFrame, snr_pred: pd.DataFrame, thresholds: pd.DataFrame, correction: str,
                      scales: Sequence[float] = tuple(np.geomspace(0.5, 4.0, 36)), scatter: float = 1.0,
                      config: Optional[str] = None, far: Optional[float] = None) -> Tuple[float, pd.DataFrame]:
    """
    ``snr_scale`` that maximises the binomial likelihood of the measured detections (``measured_efficiency``
    rows of ``correction``, one ``config``, one ``far``; defaults: the highest-k configuration present with
    k <= 2, the largest FAR) given the prediction. Returns (best scale, DataFrame scale / loglike).
    """
    m = eff_meas[eff_meas["correction"] == correction]
    far = float(m["far_per_hr"].max()) if far is None else far
    if config is None:
        config = sorted((c for c in m["config"].unique() if "of" in c), key=lambda c: int(c.split("of")[0]))[0]
    m = m[(m["config"] == config) & (m["far_per_hr"] == far)]
    sp = snr_pred[snr_pred["correction"] == correction]
    rows = []
    for sc in scales:
        pred = predicted_efficiency(sp, thresholds[thresholds["far_per_hr"] == far], {correction: (sc, scatter)})
        j = m.merge(pred[pred["config"] == config], on=["night", "observation", "correction", "radius_m", "config",
                                                          "far_per_hr"], how="inner")
        p = np.clip(j["eff_pred"].to_numpy(), 1e-4, 1.0 - 1e-4)
        ll = float(np.sum(j["hits"] * np.log(p) + (j["n"] - j["hits"]) * np.log(1.0 - p)))
        rows.append({"scale": float(sc), "loglike": ll})
    curve = pd.DataFrame(rows)
    return float(curve.loc[curve["loglike"].idxmax(), "scale"]), curve


def excess_scaling_test(profiles: Sequence[NoiseProfile], band_hz: Tuple[float, float] = SIGNAL_BAND_HZ
                        ) -> pd.DataFrame:
    """
    Does the excess noise scale with the star ("star": slope 0) or with the total signal ("total": slope -2)?
    Regresses log E on log star_fraction, E the excess PSD integrated over the signal band (star units),
    pooled with one offset per telescope and for each telescope kind alone. Rows: group, slope, slope_err, n.
    """
    pts = []
    for p in profiles:
        sel = (p.freq >= band_hz[0]) & (p.freq <= min(band_hz[1], p.freq[-1]))
        e = _trapz(p.S_excess[sel], p.freq[sel]) if sel.sum() > 1 else 0.0
        if e > 0.0 and 0.0 < p.star_fraction <= 1.0:
            pts.append((p.telescope, np.log(p.star_fraction), np.log(e)))
    out = []

    def fit(sub, label):
        n = len(sub)
        if n < 3:
            return
        tels = sorted({t for t, _, _ in sub})
        x = np.array([a for _, a, _ in sub])
        y = np.array([b for _, _, b in sub])
        cols = [x] + [np.array([t == k for t, _, _ in sub], float) for k in tels]
        A = np.column_stack(cols)
        if n <= A.shape[1] or np.ptp(x) == 0.0:
            return
        coef, res, rank, _ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ coef
        cov = np.linalg.pinv(A.T @ A) * float(resid @ resid) / (n - A.shape[1])
        out.append({"group": label, "slope": float(coef[0]), "slope_err": float(np.sqrt(cov[0, 0])), "n": n})

    fit(pts, "pooled")
    fit([q for q in pts if q[0].startswith("MAGIC")], "MAGIC")
    fit([q for q in pts if q[0].startswith("LST")], "LST")
    return pd.DataFrame(out, columns=["group", "slope", "slope_err", "n"])


def efficiency_r50(eff: pd.DataFrame, column: str, by: Sequence[str]) -> pd.DataFrame:
    """R50 (``r50``, status) of ``column`` (an efficiency vs radius_m) for every group of ``by``."""
    rows = []
    for key, g in eff.groupby(list(by)):
        g = g.sort_values("radius_m")
        value, status = r50(g["radius_m"].to_numpy(), g[column].to_numpy())
        rows.append({**dict(zip(by, key if isinstance(key, tuple) else (key,))), "r50_m": value, "status": status})
    return pd.DataFrame(rows)
