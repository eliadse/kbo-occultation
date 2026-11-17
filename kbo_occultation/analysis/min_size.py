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
import warnings
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
from ..simulation.physics import AU_m, mas_to_rad, star_radius_from_mag, teff_from_spectral_type
from .bands import SIGNAL_BAND_HZ
from .coincidence import coincidence_tolerance_s
from .spectra import _trapz, blockwise_welch, lightcurve_psd, matched_filter_sigma

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


def generic_targets(mags: Optional[Sequence[float]] = None, ref_mag: float = 7.0,
                    ref_radius_mas: float = 0.03, temperature_K: float = 10000.0) -> List[TargetStar]:
    """B = 7 ... 9 (step 0.25) stars whose radius shrinks with magnitude (star_radius_from_mag)."""
    mags = np.arange(7.0, 9.0 + 1e-9, 0.25) if mags is None else mags
    return [TargetStar(f"B{m:.2f}", float(m), star_radius_from_mag(m, ref_mag, ref_radius_mas), temperature_K)
            for m in mags]


CATALOGUE_COLUMNS = {"name": "name", "mag_B": "mag_B", "angular_radius_mas": "angular_radius_mas",
                     "temperature_K": "temperature_K", "spectral_type": "spectral_type"}


def parse_column_map(text: Optional[str]) -> Optional[Dict[str, str]]:
    """``"name=Name,mag_B=Bmag"`` -> {"name": "Name", "mag_B": "Bmag"} for ``load_target_catalogue`` (None if empty)."""
    if not text:
        return None
    out = {}
    for item in text.split(","):
        key, sep, value = item.partition("=")
        key, value = key.strip(), value.strip()
        if not sep or not value:
            raise ValueError(f"bad column mapping {item!r}: expected field=column")
        if key not in CATALOGUE_COLUMNS:
            raise ValueError(f"unknown catalogue field {key!r} (one of {', '.join(CATALOGUE_COLUMNS)})")
        out[key] = value
    return out


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
    radii_m: np.ndarray = field(default_factory=lambda: np.geomspace(10.0, 1000.0, 41))
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
    time_resolved: bool = True           # noise level per block of one search sigma-window (profiles need blocks)
    block_nperseg: int = 1024            # Welch segment of the block PSDs [binned samples]

    @property
    def bin_dt_s(self) -> float:
        return self.bin_samples * FRAME_S

    @property
    def block_s(self) -> float:
        """Length of the blocks of the time-resolved noise: the search's sigma window, 20 template durations."""
        return 20.0 * 2.0 * self.x_half_m / self.velocity_mps()

    def calibrated(self) -> "MinSizeConfig":
        """
        Copy with ``snr_scale`` / ``snr_scatter`` of this correction from ``SNR_CALIBRATION`` (time-resolved noise) or
        ``SNR_CALIBRATION_STATIONARY`` (unchanged if absent).
        """
        table = SNR_CALIBRATION if self.time_resolved else SNR_CALIBRATION_STATIONARY
        scale, scatter = table.get(self.correction, (self.snr_scale, self.snr_scatter))
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
    r_max_0 = 1.2 * np.hypot(half + r_star, r_star)
    band = BandpassConfig(*instrument.response_band())
    # the radial grid widens with the impact parameter; its step stays the b = 0 one
    n_r_grid = int(round(cfg.n_r_grid * r_max / r_max_0))
    return OccultationEngine(StarConfig(star.temperature_K, star.angular_radius_mas), band, grid,
                             NumericalConfig(n_r_grid=n_r_grid, n_star_side=cfg.n_star_side),
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


def _fine_time(cfg: MinSizeConfig) -> np.ndarray:
    """Frame-cadence time axis of the shared template span, a whole number of bins long."""
    v = cfg.velocity_mps()
    t_fine = np.arange(-cfg.x_half_m / v, cfg.x_half_m / v, FRAME_S)
    return t_fine[:len(t_fine) // cfg.bin_samples * cfg.bin_samples]


def _binned_residual(t_axis_s: np.ndarray, intensity: np.ndarray, t_fine: np.ndarray, cfg: MinSizeConfig) -> np.ndarray:
    """I - 1 of a pattern (time axis ``t_axis_s``) at the bin cadence on the span of ``t_fine``."""
    fine = np.interp(t_fine, t_axis_s, intensity)
    return fine.reshape(-1, cfg.bin_samples).mean(axis=1) - 1.0


def template_family(star: TargetStar, cfg: MinSizeConfig, instrument) -> Dict[float, np.ndarray]:
    """
    {radius_m: I - 1} at the bin cadence, for ``star`` seen through ``instrument.response()``
    (one engine, the Fresnel profiles are shared by the radii that repeat across calls).
    """
    v = cfg.velocity_mps()
    t_fine = _fine_time(cfg)
    return {radius: _binned_residual(x_m / v, intensity, t_fine, cfg)
            for radius, (x_m, intensity) in injection_patterns(star, cfg, instrument).items()}


def bank_residuals(cfg: MinSizeConfig, mimic_cfg=None, star: Optional[TargetStar] = None) -> Dict[float, np.ndarray]:
    """
    {radius_m: I - 1} of the templates of the real search (``mimic.search_templates``: the star through the
    320-450 nm band, at the shadow velocity of the search) at the bin cadence of ``template_family``. The star is the
    fixed 0.03 mas one of ``mimic_cfg`` unless ``star`` is given: the bank of that star (size and temperature), as
    ``mimic.star_bank`` builds it for an observed one.
    """
    from ..simulation.kinematics import EARTH_SEMIMAJOR_AXIS_AU
    from .mimic import MimicConfig, search_templates

    mimic_cfg = mimic_cfg or MimicConfig()
    v_bank = shadow_velocity_opposition(mimic_cfg.distance_au + EARTH_SEMIMAJOR_AXIS_AU)
    t_fine = _fine_time(cfg)
    templates = (search_templates(mimic_cfg, v_bank) if star is None else
                 search_templates(mimic_cfg, v_bank, star.angular_radius_mas, star.temperature_K))
    return {radius: _binned_residual(t, intensity, t_fine, cfg) for radius, (t, intensity) in templates.items()}


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
    block_freq: Optional[np.ndarray] = None        # Hz, grid of the block PSDs
    S_excess_blocks: Optional[np.ndarray] = None   # (n_block, len(block_freq)) as S_excess; NaN = block not usable
    white_sigma: Optional[float] = None            # star-unit sigma per bin replacing the noise model's (synthetic runs)

    @property
    def n_block(self) -> int:
        return 0 if self.S_excess_blocks is None else int(self.S_excess_blocks.shape[0])

    @property
    def fs(self) -> float:
        return 1.0 / self.dt

    _META = ("night", "observation", "channel", "telescope", "mag_B", "teff_K", "dt", "star_fraction",
             "star_fraction_model", "sf_source", "airmass", "psf_containment", "throughput", "k_B",
             "clipped_fraction", "live_s", "white_sigma")


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
                           include_excluded: bool = False, lcs=None, nperseg: int = 4096, blocks: bool = True,
                           progress: Optional[Callable[[str], None]] = None) -> List[NoiseProfile]:
    """
    One NoiseProfile per (star run, channel) of the selected nights: the run is cropped to its
    good-time window, binned (``cfg.bin_samples``), and its Welch PSD converted to star units.
    ``nperseg`` is the segment length in binned samples (4096 = 21 s, 0.05 Hz resolution).
    With ``blocks`` the run is also cut into blocks of about ``cfg.block_s`` (the search's sigma window) and
    each block gets its own median-Welch excess PSD (``S_excess_blocks``), for ``cfg.time_resolved``.
    The expensive step is loading the light curves (``lcs`` accepts the dict of
    ``snr.load_lightcurves``); use save_profiles to keep the result.
    """
    from ..observations.cc_reports import conditions_table
    from ..observations.photometry import good_time_lightcurve
    from .mimic import bin_lightcurve
    from .spectra import fractional_flux
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
                if blocks:
                    span = float(lc.time[-1] - lc.time[0]) + profile.dt
                    n_block = max(1, int(round(span / cfg.block_s)))
                    if span / n_block >= cfg.block_nperseg * profile.dt:
                        profile.block_freq, b_psd = blockwise_welch(lc.time, fractional_flux(lc.signal), profile.fs,
                                                                    n_block, cfg.block_nperseg)
                        profile.S_excess_blocks = np.maximum(b_psd / sf ** 2 - white, 0.0)
                profiles.append(profile)
                say(f"{night.name} {obs.name} {ch}: sf={sf:.3f} ({profile.sf_source}), "
                    f"clipped {profile.clipped_fraction:.2f}")
    return profiles


def load_or_measure_profiles(log, nights, cfg: Optional[MinSizeConfig] = None, path: Optional[str] = None, *,
                             include_excluded: bool = False, blocks: bool = True, force: bool = False,
                             progress: Optional[Callable[[str], None]] = None) -> List[NoiseProfile]:
    """
    The profiles of ``nights`` from ``path`` (default outputs/min_size/profiles_<nights>.npz) if it exists, else
    measured (``measure_noise_profiles``, slow) and saved there. A file without block PSDs is measured again
    when ``blocks`` is wanted; ``force`` measures again anyway.
    """
    import os

    say = progress or (lambda m: None)
    nights = list(nights) if nights else []
    path = path or os.path.join("outputs", "min_size", "profiles_" + "_".join(nights) + ".npz")
    if os.path.exists(path) and not force:
        profiles = load_profiles(path)
        if not blocks or all(p.n_block > 0 for p in profiles):
            say(f"profiles: {len(profiles)} from {path}")
            return profiles
        say(f"{path} has no block PSDs; measuring again")
    say("measuring noise profiles (slow: loads the light curves)")
    profiles = measure_noise_profiles(log, nights or None, cfg, include_excluded=include_excluded, blocks=blocks,
                                      progress=say)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    save_profiles(profiles, path)
    return profiles


def save_profiles(profiles: Sequence[NoiseProfile], path: str) -> None:
    meta = [{k: getattr(p, k) for k in NoiseProfile._META} for p in profiles]
    arrays = {}
    for i, p in enumerate(profiles):
        arrays[f"freq_{i}"], arrays[f"exc_{i}"] = p.freq, p.S_excess
        if p.S_excess_blocks is not None:
            arrays[f"bfreq_{i}"], arrays[f"bexc_{i}"] = p.block_freq, p.S_excess_blocks
    np.savez_compressed(path, meta=json.dumps(meta), **arrays)


def load_profiles(path: str) -> List[NoiseProfile]:
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z["meta"]))
        return [NoiseProfile(freq=z[f"freq_{i}"], S_excess=z[f"exc_{i}"],
                             block_freq=z[f"bfreq_{i}"] if f"bfreq_{i}" in z.files else None,
                             S_excess_blocks=z[f"bexc_{i}"] if f"bexc_{i}" in z.files else None, **m)
                for i, m in enumerate(meta)]


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
                 gain: Optional[np.ndarray] = None, signal: Optional[np.ndarray] = None) -> float:
    """
    Expected matched-filter SNR of the template (I - 1, sampled at 1 / fs) in the noise PSD ``psd``
    (one-sided, 1/Hz, on the rfft grid ``freq`` of n_fft = 2 (len(freq) - 1) points, as
    spectra.segmented_welch) after the correction ``gain`` (None = 1). For white noise of sigma s
    and no correction it is |T| / s.

    ``signal`` (I - 1 of a different event, same sampling): the template does not match it, the
    output is the largest correlation over lag, <G t, G s> / |G t|, over the same noise, so a
    signal equal to the template gives the matched value and any other of the same norm less.
    """
    n_fft = 2 * (len(freq) - 1)
    t = np.asarray(template_residual, dtype=float)
    if len(t) > n_fft:
        raise ValueError(f"template ({len(t)} samples) longer than the PSD's FFT ({n_fft})")
    g = np.ones_like(freq) if gain is None else np.asarray(gain, float)
    t_fft = g * np.fft.rfft(t, n=n_fft)
    t_corr = np.fft.irfft(t_fft, n=n_fft)
    sigma = matched_filter_sigma(freq, g ** 2 * np.asarray(psd, float), t_corr, fs)
    norm = float(np.linalg.norm(t_corr))
    if signal is None:
        return norm / sigma
    sig = np.asarray(signal, dtype=float)
    if len(sig) > n_fft:
        raise ValueError(f"signal ({len(sig)} samples) longer than the PSD's FFT ({n_fft})")
    s_fft = g * np.fft.rfft(sig, n=n_fft)
    xcorr = np.fft.irfft(np.conj(t_fft) * s_fft, n=n_fft)
    return float(np.max(xcorr) / (norm * sigma))


def bank_snr(signal: np.ndarray, bank: Mapping[float, np.ndarray], freq: np.ndarray, psd: np.ndarray, fs: float,
             gain: Optional[np.ndarray] = None) -> float:
    """Expected SNR of ``signal`` found with the best template of ``bank`` ({radius: I - 1}), as the search does."""
    return max(expected_snr(t, freq, psd, fs, gain, signal=signal) for t in bank.values())


def noise_psd(profile: NoiseProfile, star: TargetStar, cfg: MinSizeConfig, *, theory_only: bool = False,
              focus_gain: float = 1.0, block: Optional[int] = None) -> np.ndarray:
    """
    Noise PSD (star units, 1/Hz, on profile.freq) of ``star`` in the conditions of ``profile``: the
    white level of the noise model plus the run's excess (not in ``theory_only``). ``block`` selects
    the excess of that time block instead (on profile.block_freq). ``profile.white_sigma``, when set, replaces the
    noise model's sigma per bin (synthetic runs at a given SNR).
    """
    budget = model_budget(profile, star.mag_B, star.temperature_K, focus_gain)
    freq = profile.freq if block is None else profile.block_freq
    sigma = float(budget["sigma_total"]) if profile.white_sigma is None else float(profile.white_sigma)
    psd = np.full_like(freq, 2.0 * sigma ** 2 / profile.fs)
    if theory_only:
        return psd
    excess = profile.S_excess if block is None else profile.S_excess_blocks[block]
    if cfg.excess_scaling == "total":
        excess = excess * (profile.star_fraction_model / _star_fraction(budget)) ** 2
    elif cfg.excess_scaling != "star":
        raise ValueError(f"unknown excess_scaling {cfg.excess_scaling!r}")
    return psd + excess


def snr_evaluator(star: TargetStar, profile: NoiseProfile, cfg: MinSizeConfig, *, theory_only: bool = False,
                  focus_gain: float = 1.0, time_resolved: Optional[bool] = None
                  ) -> Callable[..., np.ndarray]:
    """
    ``f(template, signal=None)``: expected SNR of ``template`` (I - 1 at the bin cadence), or of ``signal``
    searched with it (``expected_snr``), in the noise of ``star`` under the conditions of ``profile``; the noise
    PSD and correction are built once. With measured noise it is multiplied by cfg.snr_scale. The result is a
    float, or (n_block,) with time-resolved noise (``per_telescope_snr``).
    """
    if not np.isclose(profile.dt, cfg.bin_dt_s, rtol=1e-3):
        raise ValueError(f"profile cadence {profile.dt} s differs from the template's {cfg.bin_dt_s} s")
    want_blocks = (cfg.time_resolved if time_resolved is None else time_resolved) and not theory_only
    use_blocks = want_blocks and profile.n_block > 0
    if want_blocks and not use_blocks:
        warnings.warn(f"{profile.night}/{profile.observation}/{profile.channel}: time-resolved noise requested but the "
                      "profile has no block PSDs (measure_noise_profiles(blocks=True)); using one PSD for the run",
                      stacklevel=2)
    psd = noise_psd(profile, star, cfg, theory_only=theory_only, focus_gain=focus_gain)
    gain = correction_gain("raw" if theory_only else cfg.correction, profile.freq, psd, profile.fs)
    scale = 1.0 if theory_only else cfg.snr_scale
    if not use_blocks:
        return lambda template, signal=None: scale * expected_snr(template, profile.freq, psd, profile.fs, gain,
                                                                  signal=signal)
    gain_b = np.interp(profile.block_freq, profile.freq, gain)
    psds_b = [None if not np.all(np.isfinite(profile.S_excess_blocks[k])) else
              noise_psd(profile, star, cfg, focus_gain=focus_gain, block=k) for k in range(profile.n_block)]

    def evaluate(template, signal=None):
        return scale * np.array([np.nan if q is None else
                                 expected_snr(template, profile.block_freq, q, profile.fs, gain_b, signal=signal)
                                 for q in psds_b])
    return evaluate


def per_telescope_snr(star: TargetStar, profile: NoiseProfile, templates: Mapping[float, np.ndarray],
                      cfg: MinSizeConfig, *, theory_only: bool = False, focus_gain: float = 1.0,
                      time_resolved: Optional[bool] = None,
                      signals: Optional[Mapping[float, np.ndarray]] = None) -> np.ndarray:
    """
    Expected SNR at every radius of cfg.radii_m. With measured noise it is multiplied by
    cfg.snr_scale (see ``SNR_CALIBRATION``: the real search normalises with a robust scatter, which on
    intermittent red noise is below the rms that the PSD gives); ``theory_only`` (white noise, no
    correction) is not scaled, the scale was fitted on real data. Stationary Gaussian noise is assumed.

    ``time_resolved`` (default ``cfg.time_resolved``; needs ``profile.S_excess_blocks``): the noise level of
    each block of the run (the search normalises with a local scatter) gives one row per block, so the
    result is (n_block, n_radii), NaN rows for unusable blocks. The correction G is still that of the whole run
    (the whitening of the search is fitted on it).

    ``signals`` {radius: I - 1}: the events that are searched (default the templates themselves, i.e. a perfect
    match); see ``expected_snr``. Used for events at an impact parameter.
    """
    f = snr_evaluator(star, profile, cfg, theory_only=theory_only, focus_gain=focus_gain, time_resolved=time_resolved)
    out = np.array([f(templates[float(r)], None if signals is None else signals[float(r)]) for r in cfg.radii_m])
    return out if out.ndim == 1 else out.T


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

def array_threshold(curves: Mapping[str, RateCurve], channels: Sequence[str], k: int, far_per_hr: float,
                    tolerance_s: float) -> Tuple[float, bool]:
    """Per-telescope SNR threshold at ``far_per_hr`` for ``k`` of the ``channels``: single_threshold if k = 1, else
    coincidence_threshold. Returns (threshold, extrapolated)."""
    sub = {ch: curves[ch] for ch in channels}
    if k == 1:
        return single_threshold(sub[channels[0]], far_per_hr)
    return coincidence_threshold(sub, tolerance_s, far_per_hr, k)


def detection_probability(snr_by_tel: Mapping[str, np.ndarray], threshold, min_ntel: int,
                          snr_scatter: float = 1.0) -> np.ndarray:
    """
    P(at least ``min_ntel`` telescopes above their threshold) per radius, each SNR Gaussian around its
    expected value. ``threshold``: a number, or {telescope: number}.

    A telescope's SNR is (n_radii,) or, with time-resolved noise, (n_block, n_radii): the probability is
    then computed per block (the same event time in every telescope) and averaged over the blocks that are
    usable in all of them, i.e. the event falls at a uniform time of the run.
    """
    arrs = {tel: np.asarray(snr, float) for tel, snr in snr_by_tel.items()}
    if any(a.ndim == 2 for a in arrs.values()):
        n_block = min(a.shape[0] if a.ndim == 2 else 1 for a in arrs.values())
        stack = np.stack([a[:n_block] if a.ndim == 2 else np.broadcast_to(a, (n_block, a.shape[-1]))
                          for a in arrs.values()])                        # (n_tel, n_block, n_radii)
        ok = np.all(np.isfinite(stack), axis=(0, 2))
        if not ok.any():
            raise ValueError("no time block is usable in every telescope")
        tels = list(arrs)
        per_block = [detection_probability({t: stack[i, b] for i, t in enumerate(tels)}, threshold, min_ntel,
                                           snr_scatter) for b in np.flatnonzero(ok)]
        return np.mean(per_block, axis=0)
    probs = []
    for tel, snr in arrs.items():
        thr = threshold[tel] if isinstance(threshold, Mapping) else threshold
        probs.append(ndtr((snr - thr) / snr_scatter))
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

def scale_quantiles(per_run_csv: str, correction: str, quantiles: Tuple[float, float] = (0.16, 0.84),
                    sf_mode: str = "measured") -> Tuple[float, float]:
    """
    (low, high) quantiles of the per-run, per-channel SNR scale of ``correction`` in the calibration's
    ``snr_calibration_per_run.csv``: the run-to-run spread of the calibration, for ``min_size_table(scale_band=)``.
    """
    df = pd.read_csv(per_run_csv)
    df = df[(df["correction"] == correction) & (df["sf_mode"] == sf_mode)].dropna(subset=["scale"])
    if df.empty:
        raise ValueError(f"{per_run_csv}: no calibration rows for {correction!r} / {sf_mode!r}")
    lo, hi = np.quantile(df["scale"].to_numpy(), quantiles)
    return float(lo), float(hi)


def min_size_table(targets: Sequence[TargetStar], profiles: Sequence[NoiseProfile], cfg: MinSizeConfig,
                   curves: Mapping[str, RateCurve], *, noise_models: Sequence[str] = ("measured", "theory"),
                   focus_gains: Sequence[float] = (1.0,), excess_scalings: Optional[Sequence[str]] = None,
                   scale_band: Optional[Tuple[float, float]] = None, return_efficiency: bool = False,
                   progress: Optional[Callable[[str], None]] = None):
    """
    R50 for every target x run (night, observation) x noise model x LST focus gain x excess scaling x
    telescope configuration (each telescope alone, ``k``of``n`` of the run's telescopes) x false-alarm rate
    (cfg.coinc_far_per_hr, then cfg.extra_fars).

    ``curves``: {channel: RateCurve} (rate_curves). The thresholds are those of the real data for both
    noise models, so ``theory`` isolates the effect of the noise on the SNR. Single telescopes use the
    same rate per telescope-hour as the coincidence FAR; arrays use coincidence_threshold.

    ``excess_scalings`` (default ``(cfg.excess_scaling,)``): the measured noise is evaluated for each of them
    ("star" and "total" bracket how the excess noise follows a change of the star light; ``theory`` has none).
    ``scale_band`` = (low, high) absolute values of ``snr_scale`` (``scale_quantiles``): ``r50_lo_m`` / ``r50_hi_m``
    are the R50 at the high / low scale, the calibration's run-to-run spread (measured noise only).

    Rows: target, mag_B, star_mas, night, observation, noise, excess_scaling ("none" for theory), lst_focus, config,
    far_per_hr, threshold, thr_extrapolated, thr_at_floor (the threshold is the candidate floor: the FAR is then
    only an upper bound, the R50 conservative), r50_m, r50_status, r50_lo_m, r50_hi_m, snr_100m (mean over the
    telescopes, interpolated; the median over the time blocks when the noise is time resolved).
    With ``return_efficiency`` also returns the long table of the detection probability per radius (the same keys
    without the r50 columns, plus radius_m and eff).
    """
    from ..simulation.instruments import Instrument

    say = progress or (lambda m: None)
    fars = [cfg.coinc_far_per_hr, *cfg.extra_fars]
    tol = cfg.tolerance_s()
    scalings = tuple(excess_scalings) if excess_scalings else (cfg.excess_scaling,)
    runs: Dict[Tuple[str, str], List[NoiseProfile]] = {}
    for p in profiles:
        runs.setdefault((p.night, p.observation), []).append(p)
    instruments = {kind: Instrument(kind, filter_type=cfg.filter_type, site="LaPalma") for kind in ("MAGIC", "LST")}

    # thresholds depend only on the telescope set and the FAR
    thr_cache: Dict[tuple, Tuple[float, bool]] = {}

    def threshold(channels: Tuple[str, ...], k: int, far: float):
        key = (channels, k, far)
        if key not in thr_cache:
            thr_cache[key] = array_threshold(curves, channels, k, far, tol)
        return thr_cache[key]

    def stretch(snr, factor):
        return {c: v * factor for c, v in snr.items()}

    rows, eff_rows = [], []
    for star in targets:
        templates = {kind: template_family(star, cfg, instr) for kind, instr in instruments.items()}
        say(f"{star.name}: B={star.mag_B:.2f}, r*={star.angular_radius_mas:.4f} mas")
        for (night, obs), group in runs.items():
            group = [p for p in group if p.channel in curves]
            for noise in noise_models:
                for scaling in (scalings if noise == "measured" else ("none",)):
                    cfg_s = cfg if scaling == "none" else replace(cfg, excess_scaling=scaling)
                    for gain in focus_gains:
                        if gain != 1.0 and not any(p.telescope.startswith("LST") for p in group):
                            continue
                        snr = {}
                        for p in group:
                            kind = "LST" if p.telescope.startswith("LST") else "MAGIC"
                            snr[p.channel] = per_telescope_snr(star, p, templates[kind], cfg_s,
                                                               theory_only=(noise == "theory"), focus_gain=gain)
                        chans = tuple(snr)
                        configs = [(p.telescope, (p.channel,), 1) for p in group]
                        configs += [(f"{k}of{len(chans)}", chans, k) for k in range(2, len(chans) + 1)]
                        typical = {c: (np.nanmedian(v, axis=0) if v.ndim == 2 else v) for c, v in snr.items()}
                        for label, sel, k in configs:
                            mean_snr = np.mean([typical[c] for c in sel], axis=0)
                            for far in fars:
                                thr, extrap = threshold(sel, k, far)
                                chosen = {c: snr[c] for c in sel}
                                prob = detection_probability(chosen, thr, k, cfg.snr_scatter)
                                value, status = r50(cfg.radii_m, prob)
                                lo = hi = float("nan")
                                if scale_band is not None and noise == "measured":
                                    lo = r50(cfg.radii_m, detection_probability(
                                        stretch(chosen, scale_band[1] / cfg.snr_scale), thr, k, cfg.snr_scatter))[0]
                                    hi = r50(cfg.radii_m, detection_probability(
                                        stretch(chosen, scale_band[0] / cfg.snr_scale), thr, k, cfg.snr_scatter))[0]
                                key = {"target": star.name, "mag_B": star.mag_B, "star_mas": star.angular_radius_mas,
                                       "night": night, "observation": obs, "noise": noise,
                                       "excess_scaling": scaling, "lst_focus": gain, "config": label,
                                       "far_per_hr": far}
                                rows.append({**key, "threshold": thr, "thr_extrapolated": extrap,
                                             "thr_at_floor": bool(thr <= SNR_FLOOR + 1e-9),
                                             "r50_m": value, "r50_status": status, "r50_lo_m": lo, "r50_hi_m": hi,
                                             "snr_100m": float(np.interp(np.log(100.0), np.log(cfg.radii_m),
                                                                         mean_snr))})
                                if return_efficiency:
                                    eff_rows.extend({**key, "radius_m": float(r), "eff": float(e)}
                                                    for r, e in zip(cfg.radii_m, prob))
    table = pd.DataFrame(rows)
    return (table, pd.DataFrame(eff_rows)) if return_efficiency else table


def summarize(table: pd.DataFrame) -> pd.DataFrame:
    """Median / min / max R50 over the runs, per target, noise, excess scaling, focus, configuration and FAR."""
    keys = ["target", "mag_B", "star_mas", "noise", "excess_scaling", "lst_focus", "config", "far_per_hr"]
    keys = [k for k in keys if k in table]
    agg = dict(r50_median=("r50_m", "median"), r50_min=("r50_m", "min"), r50_max=("r50_m", "max"),
               n_runs=("r50_m", "size"), n_unreached=("r50_status", lambda s: int((s == "above_grid").sum())),
               n_bound=("r50_status", lambda s: int((s == "below_grid").sum())))
    if "r50_lo_m" in table:
        agg.update(r50_lo_median=("r50_lo_m", "median"), r50_hi_median=("r50_hi_m", "median"))
    if "thr_at_floor" in table:
        agg.update(thr_at_floor=("thr_at_floor", "max"))
    return table.groupby(keys).agg(**agg).reset_index()


# ───────────────────────────── impact parameter ─────────────────────────────

# impact parameters are sampled at these fractions of H(R) / 2; the last nodes catch the diffraction fringes
# that stay detectable outside the geometric shadow of a bright star
B_FRACTIONS = (0.0, 0.2, 0.4, 0.6, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 1.8, 2.2, 2.8, 3.5, 4.5, 6.0)
B_MAX_M = 15000.0
WIDTH_WAVELENGTH_M = 420e-9           # for H(r) (Fresnel scale), as tests/observing_time_vs_kbo_size.py


def shadow_width_m(radius_m, star: TargetStar, cfg: MinSizeConfig, wavelength_m: float = WIDTH_WAVELENGTH_M):
    """Effective shadow width H(r) of Nihei et al. 2007 (``upper_limits.effective_shadow_width_m``) of ``star``."""
    from .upper_limits import effective_shadow_width_m

    r_star = star.angular_radius_mas * mas_to_rad * cfg.distance_au * AU_m
    return effective_shadow_width_m(radius_m, cfg.distance_au, wavelength_m, r_star)


def impact_nodes(radius_m: float, star: TargetStar, cfg: MinSizeConfig, wavelength_m: float = WIDTH_WAVELENGTH_M
                 ) -> np.ndarray:
    """Impact parameters [m] at which the detection probability of radius ``radius_m`` is sampled."""
    nodes = 0.5 * float(shadow_width_m(radius_m, star, cfg, wavelength_m)) * np.asarray(B_FRACTIONS)
    return nodes[nodes <= B_MAX_M]


def width_signals(star: TargetStar, cfg: MinSizeConfig, instrument, radii: Optional[Sequence[float]] = None,
                  wavelength_m: float = WIDTH_WAVELENGTH_M) -> Dict[float, Tuple[np.ndarray, np.ndarray]]:
    """
    {radius_m: (b_m (n_b,), I - 1 (n_b, n_bin))}: the event of ``star`` seen through ``instrument`` at the impact
    parameters ``impact_nodes``, on the bin cadence of ``template_family``. One engine (its radial grid covers the
    largest impact parameter) and the Fresnel profile of each radius is computed once.
    """
    radii = [float(r) for r in (cfg.radii_m if radii is None else radii)]
    nodes = {r: impact_nodes(r, star, cfg, wavelength_m) for r in radii}
    b_top = max(float(n.max()) for n in nodes.values())
    wide = _engine(star, replace(cfg, impact_parameter_m=b_top), instrument)
    v = cfg.velocity_mps()
    t_fine = _fine_time(cfg)
    out = {}
    for r in radii:
        rows = []
        for b in nodes[r]:
            x_m, intensity = wide.compute(KBOConfig(r, cfg.distance_au, float(b)), star.angular_radius_mas)
            rows.append(_binned_residual(x_m / v, intensity, t_fine, cfg))
        out[r] = (nodes[r], np.array(rows))
    return out


def detection_width(star: TargetStar, group: Sequence[NoiseProfile], cfg: MinSizeConfig,
                    templates: Mapping[str, Mapping[float, np.ndarray]],
                    signals: Mapping[str, Mapping[float, Tuple[np.ndarray, np.ndarray]]], *, k: int, threshold: float,
                    focus_gain: float = 1.0, theory_only: bool = False,
                    wavelength_m: float = WIDTH_WAVELENGTH_M) -> pd.DataFrame:
    """
    Detection width W(R) = 2 integral_0^inf P(R, b) db of the telescopes ``group`` (one run, the channels that
    take part) detecting an event with ``k`` of them above ``threshold``: P(R, b) is ``detection_probability`` of
    the expected SNR of the event at impact parameter b searched with the b = 0 template (what the search does).
    ``templates`` and ``signals``: {"MAGIC" | "LST": template_family / width_signals} of ``star``, on cfg.radii_m.
    The telescope offsets (~100 m) are neglected against the impact parameters (km).

    Rows: radius_m, h_m (H(R)), width_m, eff_eff (W / H, it can exceed 1), p_b0 (P at b = 0, the efficiency used
    elsewhere), p_tail (P at the outermost node; > ~0.02 means the integral is cut too early).
    """
    evals = {p.channel: snr_evaluator(star, p, cfg, theory_only=theory_only, focus_gain=focus_gain) for p in group}
    kind = {p.channel: "LST" if _is_lst(p) else "MAGIC" for p in group}
    rows = []
    for r in cfg.radii_m:
        r = float(r)
        snr = {}
        for ch, f in evals.items():
            b, resid = signals[kind[ch]][r]
            vals = np.array([f(templates[kind[ch]][r], resid[i]) for i in range(len(b))])    # (n_b,) or (n_b, n_block)
            snr[ch] = vals if vals.ndim == 1 else vals.T
        prob = detection_probability(snr, threshold, k, cfg.snr_scatter)
        width = 2.0 * _trapz(prob, b)
        h = float(shadow_width_m(r, star, cfg, wavelength_m))
        rows.append({"radius_m": r, "h_m": h, "width_m": width, "eff_eff": width / h, "p_b0": float(prob[0]),
                     "p_tail": float(prob[-1])})
    return pd.DataFrame(rows)


# ───────────────────────────── calibration against injections ───────────────

# correction -> (snr_scale, snr_scatter) fitted on array injections (examples/studies/min_size_calibration.py,
# outputs/min_size/calibration/snr_calibration.csv). Median over runs of (measured peak SNR) / (prediction); scatter
# is the absolute 1.4826 MAD of the peak SNR below SNR 20.
# SNR_CALIBRATION is for the time-resolved noise (MinSizeConfig.time_resolved, the default). 2026-11-16: 2025-12-16
# (5 single-channel runs) + 2026-03-08 + 2026-08-19 (7 runs x 3 channels), 26 run-channels, 50 trials x 12 radii
# (50-400 m), measured star fraction. Per-run-channel scale range 0.74-1.37 (whitened), 0.71-1.13 (butter_1),
# 0.37-1.23 (raw); the likelihood fit on the detections gives 1.15 / 1.02 / 0.63 (whitened / butter_1 / raw).
SNR_CALIBRATION: Dict[str, Tuple[float, float]] = {"whitened": (1.072, 1.174), "raw": (0.655, 1.774),
                                                   "butter_1": (0.971, 1.221)}
# SNR_CALIBRATION_STATIONARY is for one PSD per run (time_resolved=False; step 2, 2026-11-12: 2026-03-08 +
# 2026-08-19, 21 run-channels, per-run range 0.67-1.65 (whitened), 0.89-1.58 (butter_1), 0.87-6.5 (raw)). Its spread
# is 1.5-6 times larger than the time-resolved one, and HD170263 is 38-42 % too pessimistic with it.
SNR_CALIBRATION_STATIONARY: Dict[str, Tuple[float, float]] = {"whitened": (0.912, 1.234), "raw": (1.539, 2.250),
                                                              "butter_1": (1.082, 1.463)}

CALIBRATION_RADII_M = (50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 400.0)
_RUN = ["night", "observation", "correction", "sf_mode"]


def observed_target(obs, star_mas: float = 0.03) -> TargetStar:
    """
    The star of a logged run as a TargetStar of angular radius ``star_mas`` (a synthetic run carries its own
    ``angular_radius_mas``, which wins). ``calibrate_run`` and ``predicted_snr`` pass the star's own size,
    ``mimic.observed_star_mas`` (from its B magnitude), as the per-star search bank uses.
    """
    teff = teff_from_spectral_type(obs.spectral_type) if obs.spectral_type else None
    radius = getattr(obs, "angular_radius_mas", None) or star_mas
    return TargetStar(obs.name, float(obs.mag_B), float(radius), float(teff or 10000.0), obs.spectral_type)


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
                  patterns: Optional[dict] = None,
                  full_candidates: Optional[list] = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Array injections of one run and one correction (``noise_removal.injection_efficiency_array``).

    The search is the one that writes real_candidates.csv: the method's reshaped bank of the star's own size,
    ``mimic.star_bank(ccfg.mimic, obs, ...)`` (``bank`` replaces it, e.g. the fixed 0.03 mas bank of the older
    search). The injected pattern is the observed star, at its own size (``mimic.observed_star_mas``), through the
    telescope's own response (``injection_patterns``), diluted by the star fraction of each channel:
    ``sf_mode="measured"`` is the profile's (what the prediction uses), ``"model"`` the noise
    model's. ``patterns`` is a {kind: {R: (x, I)}} cache that is filled in (pass the same dict for the
    corrections of one run and the same ``radii``). ``full_candidates`` (a list) is extended with the candidates of
    the blind search of each uninjected channel (run keys, channel, telescope, search row): the rate curves of a
    synthetic run, which has no real_candidates.csv. Returns (trials, candidates, run_calib);
    run_calib has the full-run search statistics per channel (SNR-trace std and MAD-sigma, number of
    candidates) to compare with real_candidates.csv.
    """
    from ..observations.photometry import good_time_lightcurve
    from ..simulation.instruments import Instrument
    from ..simulation.kinematics import EARTH_SEMIMAJOR_AXIS_AU
    from .mimic import _template_duration, bin_lightcurve, observed_star_mas, search_lightcurve, star_bank
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
    templates_raw = bank if bank is not None else star_bank(ccfg.mimic, obs, v_bank)
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
        if full_candidates is not None:
            full_candidates.extend({**key, "channel": ch, "telescope": p.telescope, **r} for r in rows)
        calib_rows.append({**key, "channel": ch, "telescope": p.telescope, "n_candidates": len(rows),
                           "n_ge5": sum(r["snr"] >= 5.0 for r in rows), "snr_std": c["snr_std"],
                           "snr_mad_sigma": c["snr_mad_sigma"], "live_s": len(real.time) * real.dt})

    star = observed_target(obs, observed_star_mas(obs))                  # the star's own size, whatever the bank
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
                         "thr_at_floor": bool(thr <= SNR_FLOOR + 1e-9),
                         "timeslide_threshold": np.nan, "timeslide_resolved": np.nan})
        for k in range(2, n + 1):
            thr, ext = coincidence_threshold(curves, tolerance_s, far, k)
            ts, resolved = np.nan, np.nan
            if k == 2 and far_check is not None and len(far_check):
                ok = far_check[far_check["slide_rate_per_hr"] <= far]
                ts = float(ok["snr"].min()) if len(ok) else np.nan
                resolved = bool(far * float(far_check["slide_hours"].iloc[0]) >= 1.0)
            rows.append({"far_per_hr": far, "config": f"{k}of{n}", "k": k, "threshold": thr,
                         "thr_extrapolated": ext, "thr_at_floor": bool(thr <= SNR_FLOOR + 1e-9),
                         "timeslide_threshold": ts, "timeslide_resolved": resolved})
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
                  radii: Sequence[float] = CALIBRATION_RADII_M, bank: Optional[Mapping[float, np.ndarray]] = None
                  ) -> pd.DataFrame:
    """
    Step-1 expected SNR (``per_telescope_snr``, snr_scale 1) of the observed star of every run at the
    injection radii, per correction and channel: night, observation, correction, channel, telescope,
    radius_m, snr_pred, block. ``runs``: ``calibration_runs``. ``block`` is 0 for a stationary prediction;
    with ``cfg.time_resolved`` (and block PSDs in the profiles) every time block is a row (NaN where unusable).

    ``bank`` ({radius: I - 1} at the bin cadence, e.g. ``bank_residuals``): the event of each radius is searched
    with the best template of the bank (``bank_snr``), as the real search does, instead of its own (perfectly
    matched) template.
    """
    from ..simulation.instruments import Instrument
    from .mimic import observed_star_mas

    cfg_r = replace(cfg, radii_m=np.asarray(radii, float), snr_scale=1.0)
    instruments = {k: Instrument(k, filter_type=cfg.filter_type, site="LaPalma") for k in ("MAGIC", "LST")}
    rows = []
    for night, obs, group in runs:
        star = observed_target(obs, observed_star_mas(obs))
        kinds = {"LST" if _is_lst(p) else "MAGIC" for p in group}
        templates = {k: template_family(star, cfg_r, instruments[k]) for k in kinds}
        for corr in corrections:
            cfg_c = replace(cfg_r, correction=corr)
            for p in group:
                family = templates["LST" if _is_lst(p) else "MAGIC"]
                if bank is None:
                    snr = np.atleast_2d(per_telescope_snr(star, p, family, cfg_c))
                else:
                    f = snr_evaluator(star, p, cfg_c)
                    best = np.array([np.max([f(t, signal=family[float(r)]) for t in bank.values()], axis=0)
                                     for r in radii])                          # (n_radii,) or (n_radii, n_block)
                    snr = np.atleast_2d(best) if best.ndim == 1 else best.T
                rows.extend({"night": night.name, "observation": obs.name, "correction": corr,
                             "channel": p.channel, "telescope": p.telescope, "radius_m": float(r),
                             "snr_pred": float(v), "block": b}
                            for b, row in enumerate(snr) for r, v in zip(radii, row))
    return pd.DataFrame(rows)


def predicted_efficiency(snr_pred: pd.DataFrame, thresholds: pd.DataFrame,
                         scales: Optional[Mapping[str, Tuple[float, float]]] = None) -> pd.DataFrame:
    """
    Detection probability from ``predicted_snr`` at the same thresholds as ``measured_efficiency``, with the
    expected SNRs multiplied by ``scales[correction][0]`` and a scatter ``scales[correction][1]`` (default
    1, 1). With several ``block`` rows the probability is averaged over the blocks (``detection_probability``).
    Rows: night, observation, correction, radius_m, config, far_per_hr, threshold, eff_pred.
    """
    scales = scales or {}
    if "block" not in snr_pred:
        snr_pred = snr_pred.assign(block=0)
    thr = {(r.far_per_hr, r.config): r.threshold for r in thresholds.itertuples()}
    fars = sorted({f for f, _ in thr})
    rows = []
    for (night, obs, corr), g in snr_pred.groupby(["night", "observation", "correction"]):
        scale, scatter = scales.get(corr, (1.0, 1.0))
        radii = np.sort(g["radius_m"].unique())
        snr = {}
        for ch, gg in g.groupby("channel"):
            v = gg.pivot(index="block", columns="radius_m", values="snr_pred").sort_index()[radii].to_numpy() * scale
            snr[ch] = v[0] if len(v) == 1 else v
        n = len(snr)
        for far in fars:
            for ch, v in snr.items():
                if (far, ch) in thr:
                    p = ndtr((v - thr[(far, ch)]) / scatter)
                    p = p if p.ndim == 1 else np.nanmean(p, axis=0)
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
    if "block" in snr_pred:                                   # time-resolved prediction: the median block
        snr_pred = (snr_pred.groupby([c for c in snr_pred.columns if c not in ("block", "snr_pred")],
                                     as_index=False)["snr_pred"].median())
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
    k <= 2, else the first channel; the largest FAR) given the prediction. Returns (best scale, DataFrame scale / loglike).
    """
    m = eff_meas[eff_meas["correction"] == correction]
    far = float(m["far_per_hr"].max()) if far is None else far
    if config is None:
        multi = sorted((c for c in m["config"].unique() if "of" in c), key=lambda c: int(c.split("of")[0]))
        config = multi[0] if multi else sorted(m["config"].unique())[0]          # single-channel data: one channel
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
