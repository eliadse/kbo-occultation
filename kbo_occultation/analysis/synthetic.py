# kbo_occultation/analysis/synthetic.py
"""
Synthetic star runs for the minimum-size pipeline: Gaussian noise of a chosen level and colour, in place of a
logged run, through the same injection / blind search / threshold / efficiency / R50 code as the real calibration
(``min_size.calibrate_run``, examples/studies/synthetic_min_size.py).

Noise models (``synthetic_profile``), all in star units (fractional noise of the star signal, star level 1):

- ``white``: the white level of ``Instrument.noise_budget`` of the star in the conditions of a reference run
  (airmass, throughput, PSF), or a given SNR (``snr`` per ``snr_dt_s``);
- ``coloured``: the same white level plus the measured excess PSD of the reference run, i.e. exactly the PSD that
  ``min_size.noise_psd`` assumes for the prediction.

The noise is stationary, Gaussian and independent between channels, so a prediction that is right for these
assumptions gives an SNR scale of 1 against it; whatever differs from 1 comes from the search itself (template bank,
detrending, empirical normalisation) rather than from the real noise.

``SyntheticNight`` / ``SyntheticObservation`` have the attributes of the observation log's entries that
``calibrate_run`` and ``predicted_snr`` read, so they can stand in for real runs.
"""

import zlib
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import fft as sfft

from ..config import ORM_ARRAY, TelescopeConfig
from ..observations.photometry import LightCurve
from .bands import SIGNAL_BAND_HZ
from .min_size import (FRAME_S, MinSizeConfig, NoiseProfile, RateCurve, TargetStar, noise_psd, star_radius_from_mag)
from .spectra import _trapz

NPERSEG = 4096              # Welch segment of the profiles [binned samples], as measure_noise_profiles
T0_UNIX = 1.7e9             # start time of the synthetic runs, as a logged run's unix seconds


# ───────────────────────────── noise ────────────────────────────────────────

def gaussian_noise(n: int, fs: float, white_psd: float, freq: Optional[np.ndarray] = None,
                   excess: Optional[np.ndarray] = None, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    ``n`` samples of zero-mean Gaussian noise at sample rate ``fs`` [Hz] with the one-sided PSD [1/Hz]
    ``white_psd + excess(freq)`` (``excess`` interpolated, zero above ``freq[-1]``; no power at f = 0).
    White noise of sigma s at rate fs has ``white_psd = 2 s^2 / fs``, and the level is the same after binning,
    so the noise can be generated at the frame cadence for a PSD given at the bin cadence.
    """
    rng = rng or np.random.default_rng()
    m = sfft.next_fast_len(n, real=True)
    f = sfft.rfftfreq(m, d=1.0 / fs)
    psd = np.full(len(f), float(white_psd))
    if excess is not None:
        psd += np.interp(f, np.asarray(freq, float), np.asarray(excess, float), right=0.0)
    psd[0] = 0.0
    x = sfft.irfft(sfft.rfft(rng.standard_normal(m)) * np.sqrt(psd * fs / 2.0), n=m)
    return x[:n]


# ───────────────────────────── profiles ─────────────────────────────────────

def _band_excess(profile: NoiseProfile, band_hz: Tuple[float, float] = SIGNAL_BAND_HZ) -> float:
    sel = (profile.freq >= band_hz[0]) & (profile.freq <= min(band_hz[1], profile.freq[-1]))
    return _trapz(profile.S_excess[sel], profile.freq[sel]) if sel.sum() > 1 else 0.0


def reference_runs(profiles: Sequence[NoiseProfile], mode: str = "median", channels: Sequence[str] = ("A", "B", "C")
                   ) -> List[Dict[str, NoiseProfile]]:
    """
    Real runs whose conditions (and, for coloured noise, excess PSD) the synthetic runs copy: [{channel: profile}] of
    the runs that have every channel in ``channels``. ``mode``: ``median`` (the run with the median in-band excess,
    summed over channels), ``all`` (every complete run, by in-band excess) or ``NIGHT/OBSERVATION``.
    """
    by_run: Dict[Tuple[str, str], Dict[str, NoiseProfile]] = {}
    for p in profiles:
        by_run.setdefault((p.night, p.observation), {})[p.channel] = p
    complete = {k: g for k, g in by_run.items() if all(ch in g for ch in channels)}
    if not complete:
        raise ValueError(f"no run with all of the channels {list(channels)} among the noise profiles")
    ranked = sorted(complete.items(), key=lambda kv: sum(_band_excess(kv[1][ch]) for ch in channels))
    if mode == "median":
        return [ranked[len(ranked) // 2][1]]
    if mode == "all":
        return [g for _, g in ranked]
    key = tuple(mode.split("/", 1))
    if key not in complete:
        raise ValueError(f"reference run {mode!r} not among {sorted('/'.join(k) for k in complete)}")
    return [complete[key]]


def synthetic_profile(night: str, observation: str, channel: str, mag_B: float, cfg: MinSizeConfig, *,
                      noise: str = "white", reference: Optional[NoiseProfile] = None, snr: Optional[float] = None,
                      snr_dt_s: Optional[float] = None) -> NoiseProfile:
    """
    NoiseProfile of a synthetic run's channel. Conditions (airmass, throughput, PSF containment, k_B) are those of
    ``reference`` (a real profile), else the defaults (``cfg.airmass``, 1, 1, none). ``coloured`` takes the
    reference's excess PSD (required), ``white`` none. ``snr``: the photometric SNR of the star in bins of
    ``snr_dt_s`` seconds, which sets the white level (``NoiseProfile.white_sigma``) instead of the noise model.
    The star fraction is 1: the noise is in star units, and the injection is not diluted.
    """
    if noise not in ("white", "coloured"):
        raise ValueError(f"unknown noise model {noise!r} (white | coloured)")
    if noise == "coloured" and reference is None:
        raise ValueError("coloured noise needs a reference profile (excess PSD)")
    if snr is not None and (snr <= 0 or not snr_dt_s or snr_dt_s <= 0):
        raise ValueError("snr needs snr > 0 and snr_dt_s > 0")
    tel: TelescopeConfig = ORM_ARRAY.by_channel(channel)
    if reference is not None and reference.channel != channel:
        raise ValueError(f"reference profile is channel {reference.channel!r}, not {channel!r}")
    fs_bin = 1.0 / cfg.bin_dt_s
    if noise == "coloured":
        freq, excess = reference.freq, reference.S_excess
    else:
        freq = np.linspace(0.0, 0.5 * fs_bin, NPERSEG // 2 + 1)
        excess = np.zeros_like(freq)
    return NoiseProfile(
        night=night, observation=observation, channel=channel, telescope=tel.name, mag_B=float(mag_B),
        teff_K=10000.0, dt=cfg.bin_dt_s, freq=freq, S_excess=excess, star_fraction=1.0, star_fraction_model=1.0,
        sf_source="model", airmass=reference.airmass if reference else cfg.airmass,
        psf_containment=reference.psf_containment if reference else 1.0,
        throughput=reference.throughput if reference else 1.0, k_B=reference.k_B if reference else None,
        clipped_fraction=0.0, live_s=0.0,
        white_sigma=None if snr is None else float(np.sqrt(snr_dt_s / cfg.bin_dt_s) / snr))


# ───────────────────────────── stand-ins for the log's runs ─────────────────

@dataclass
class SyntheticNight:
    """What ``calibrate_run`` / ``predicted_snr`` read of an observation-log night."""
    name: str
    channels: List[str]
    dc_report: Optional[str] = None
    optical_filter: Optional[str] = None

    def telescope(self, channel: str) -> TelescopeConfig:
        return ORM_ARRAY.by_channel(channel)


class SyntheticObservation:
    """
    A star run of generated noise: ``lightcurve(channel)`` is ``1 + noise`` (star units) at the frame cadence,
    drawn once per channel from the PSD of ``min_size.noise_psd`` of the channel's profile (deterministic in
    ``seed`` and the run / channel names, independent between channels). ``angular_radius_mas`` is read by
    ``min_size.observed_target``.
    """
    has_cache = True
    spectral_type = None
    ra_deg = None
    dec_deg = None
    good_time = None

    def __init__(self, name: str, star: TargetStar, profiles: Mapping[str, NoiseProfile], cfg: MinSizeConfig,
                 night_name: str = "synthetic", run_s: float = 900.0, seed: int = 0, t0: float = T0_UNIX):
        self.name = name
        self.star = star
        self.mag_B = star.mag_B
        self.angular_radius_mas = star.angular_radius_mas
        self.profiles = dict(profiles)
        self.channels = sorted(self.profiles)
        self.cfg = cfg
        self.night_name = night_name
        self.run_s = float(run_s)
        self.seed = seed
        self.t0 = t0
        self.n = int(self.run_s / FRAME_S)
        self._time = self.t0 + np.arange(self.n) * FRAME_S
        self._signal: Dict[str, np.ndarray] = {}

    def time_span(self) -> Tuple[float, float]:
        return float(self._time[0]), float(self._time[-1])

    def _draw(self, channel: str) -> np.ndarray:
        p = self.profiles[channel]
        total = noise_psd(p, self.star, self.cfg)                                  # white + excess, on p.freq
        white = float(noise_psd(p, self.star, self.cfg, theory_only=True)[0])
        rng = np.random.default_rng([self.seed, zlib.crc32(f"{self.night_name}/{self.name}/{channel}".encode())])
        return 1.0 + gaussian_noise(self.n, 1.0 / FRAME_S, white, p.freq, total - white, rng)

    def lightcurve(self, channel: str, **_) -> LightCurve:
        if channel not in self._signal:
            self._signal[channel] = self._draw(channel)
        return LightCurve(self._time, self._signal[channel].copy(), meta={"channel": channel, "synthetic": True})

    def release(self) -> None:
        """Drop the generated noise (it is redrawn identically on demand)."""
        self._signal.clear()


def synthetic_runs(mags: Sequence[float], cfg: MinSizeConfig, *, noise: str = "white",
                   references: Optional[Sequence[Mapping[str, NoiseProfile]]] = None, snr: Optional[float] = None,
                   snr_dt_s: Optional[float] = None, n_runs: int = 2, run_s: float = 900.0, seed: int = 0,
                   channels: Sequence[str] = ("A", "B", "C"), ref_mag: float = 7.0, ref_radius_mas: float = 0.03):
    """
    ``n_runs`` synthetic runs per magnitude of ``mags``, as ``[(night, obs, [NoiseProfile per channel])]``
    (``min_size.calibration_runs``). The night is named ``synthetic-<noise>`` (``-snr<snr>`` with a given SNR) so
    that the run keys of the injection tables carry the noise model; the runs are ``B<mag>_<i>``. The star's
    radius follows ``min_size.star_radius_from_mag`` (as ``generic_targets``). Run ``i`` copies the conditions
    (and excess) of ``references[i % len(references)]`` ({channel: profile}, see ``reference_runs``).
    """
    if noise == "coloured" and not references:
        raise ValueError("coloured noise needs references (reference_runs)")
    night = SyntheticNight(f"synthetic-{noise}" + ("" if snr is None else f"-snr{snr:g}"), list(channels))
    out = []
    for mag in mags:
        star_radius = star_radius_from_mag(float(mag), ref_mag, ref_radius_mas)
        for i in range(n_runs):
            name = f"B{float(mag):g}_{i}"
            ref = references[i % len(references)] if references else {}
            group = [synthetic_profile(night.name, name, ch, mag, cfg, noise=noise, reference=ref.get(ch), snr=snr,
                                       snr_dt_s=snr_dt_s) for ch in channels]
            star = TargetStar(name, float(mag), star_radius)
            out.append((night, SyntheticObservation(name, star, {p.channel: p for p in group}, cfg, night.name,
                                                    run_s, seed), sorted(group, key=lambda p: p.channel)))
    return out


# ───────────────────────────── thresholds from the synthetic search ─────────

def rate_curves_from_candidates(cands: pd.DataFrame, run_calib: pd.DataFrame, n_tail: int = 20,
                                by: Sequence[str] = ("night", "correction")) -> Dict[tuple, Dict[str, RateCurve]]:
    """
    {key: {channel: RateCurve}} from the blind searches of the uninjected synthetic runs (``calibrate_run``'s
    ``full_candidates`` as a DataFrame, and its ``run_calib`` for the live time): the synthetic counterpart of
    ``min_size.rate_curves`` on real_candidates.csv. ``key`` is the tuple of the ``by`` columns, which both
    tables must have (add e.g. a magnitude column to threshold per magnitude); the live time of a channel is
    the sum of ``live_s`` over its runs (``sf_mode`` "measured" rows, as the search does not depend on it). A
    ``by`` column that is not a run key, e.g. ``radius_m`` of the matched search (one blind search per template),
    keeps one row per value of it.
    """
    by = list(by)
    calib = run_calib[run_calib["sf_mode"] == "measured"] if "sf_mode" in run_calib else run_calib
    run_keys = ["night", "observation", "correction", "channel"]
    calib = calib.drop_duplicates(run_keys + [c for c in by if c not in run_keys])
    hours = calib.groupby(by + ["channel"])["live_s"].sum() / 3600.0
    out: Dict[tuple, Dict[str, RateCurve]] = {}
    for (*key, ch), h in hours.items():
        sel = cands
        for col, val in zip(by, key):
            sel = sel[sel[col] == val] if len(sel) else sel
        sel = sel[sel["channel"] == ch] if len(sel) else sel
        snr = sel["snr"].to_numpy() if len(sel) else np.array([])
        out.setdefault(tuple(key), {})[ch] = RateCurve(snr, float(h), n_tail)
    return out
