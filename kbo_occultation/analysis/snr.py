"""
Per-night SNR report: measured SNR of every star observation in the log, by
every estimator of snr_measurement, next to the theoretical SNR of
Instrument at the observation's airmass and mirror throughput.

This is what ``kbo-night snr`` runs; in Python::

    from kbo_occultation.observations.log import load_observation_log
    from kbo_occultation.analysis.snr import snr_report

    log = load_observation_log()
    results = snr_report(log, ["2026-08-19"])   # index: night, observation, telescope, n_bin

The steps are also available one by one (``load_lightcurves``,
``nsb_table``, ``measured_snr_table``, ``theoretical_snr_table``), e.g. to
reuse the light curves for diagnostics. Nights / observations marked
``exclude`` in the log are skipped unless named or ``include_excluded=True``.
Measured-SNR columns (all at the
binned cadence, on the good-time window, NSB subtracted):

- mean_fast / peak_fast: mean / last histogram peak of the fast signal over
  the mean within-bin standard error (old notebook method 1 / p1);
- mean_comb / peak_comb: same on the fast+DC combined signal (methods 2/3);
- mean_fast_dcnsb: mean_fast with the NSB from the DC background pixel;
- detr: baseline over the robust bin-to-bin scatter after removing a
  DETREND_WINDOW_S rolling baseline -- **the one to compare with theory**;
- inband_rms / mf_sigma / welch_short: noise in the diffraction-signal band, see
  measured_snr_table;
- detr_full / detr_despiked / detr_dcdetrend: detr on the whole run, after
  dropping DC excursions, after dividing out the smoothed DC trend.

Theory uses the Johnson B magnitude with the CALSPEC Vega zero point, times a colour term from the
star's T_eff (``teff_K``, from the log's spectral type; Vega colours if it has none).
Theory columns, per band: theo_shot_<band> (shot noise only), theo_zenith_<band>
(shot x F^2 + electronic + scintillation at the zenith), theo_obs_<band> (same at
the observation's airmass with extinction, times the mirror area left after
the faulty panels).
"""

import os
from typing import Dict, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ..simulation.physics import teff_from_spectral_type
from ..observations.cc_reports import conditions_table
from ..observations.dc_combine import (combine_lightcurve_with_dc, despike_lightcurve_with_dc,
                                       detrend_lightcurve_with_dc)
from ..observations.photometry import LightCurve, good_time_lightcurve
from .bands import SIGNAL_BAND_HZ, default_template_residual
from .spectra import noise_band_numbers
from .snr_measurement import (dc_to_fast_scale, detrended_snr, measured_snr,
                              nsb_from_dark, nsb_from_dc_background, nsb_from_dc_histogram,
                              throughput_from_faulty_panels)

N_BINS = (10, 20, 200)            # samples per bin: ~5 ms, ~10ms (minimum pattern) and ~100 ms (old SNR plots)
from ..config import DC_CUT_HZ   # fast/DC crossover (DC cadence: MAGIC 1 s, LST ~2 s)
DETREND_WINDOW_S = 10.0       # rolling-median baseline for detr (~20x an occultation template)
BANDS = ("B",)

LightCurves = Dict[Tuple[str, str, str], LightCurve]   # (night, observation, channel) -> LightCurve


def load_lightcurves(log, nights=None, *, include_excluded=False) -> LightCurves:
    """Outlier-cleaned light curve of every cached observation and channel of the selected nights."""
    return {(night.name, obs.name, ch): obs.lightcurve(ch)
            for night in log.select(nights, include_excluded) for obs in night.observations.values() if obs.has_cache
            for ch in night.channels}


def _lcs(log, nights, lcs, include_excluded=False):
    return lcs if lcs is not None else load_lightcurves(log, nights, include_excluded=include_excluded)


def _has_dc(night):
    return night.dc_report is not None and os.path.exists(night.dc_report)


def nsb_table(log, nights=None, n_bin: int = N_BINS[0], lcs: Optional[LightCurves] = None, *,
              include_excluded=False) -> pd.DataFrame:
    """
    NSB (+ electronic pedestal) level per star observation and telescope
    (index night, observation, telescope), in fast (variance) units, by three
    estimators:

    - nsb_dark: mean of the night's dark run, binned by ``n_bin`` (error: nsb_dark_err);
    - nsb_dc_bkg: lower DC background pixel during the run, times dc_scale;
    - nsb_dc_hist: first peak of the on-target DC histogram, times dc_scale
      (meant for LST-1; wrong if the camera was closed during the report);

    and ``nsb_default`` / ``nsb_source``: dark if the night has one, else
    dc_hist for LST-1 and dc_bkg for MAGIC. ``star_mean`` is the mean fast
    signal in the good-time window.
    """
    lcs = _lcs(log, nights, lcs, include_excluded)
    rows = []
    for night in log.select(nights, include_excluded):
        for ch in night.channels:
            dc_key = night.telescope(ch).dc_key
            dark = night.dark if night.dark is not None and night.dark.has_cache else None
            dark_nsb = (nsb_from_dark(good_time_lightcurve(lcs[(night.name, dark.name, ch)], dark.good_time), n_bin)
                        if dark else None)
            for obs in night.stars.values():
                if not obs.has_cache:
                    continue
                lc = lcs[(night.name, obs.name, ch)]
                row = {"night": night.name, "observation": obs.name, "telescope": night.telescope(ch).name,
                       "dc_key": dc_key, "nsb_dark": dark_nsb[0] if dark_nsb else np.nan,
                       "nsb_dark_err": dark_nsb[1] if dark_nsb else 0.0, "dc_scale": np.nan}
                if _has_dc(night):
                    try:
                        scale = dc_to_fast_scale(lc, night.dc_report, dc_key)
                        bkg, _ = nsb_from_dc_background(night.dc_report, dc_key, lc.time[0], lc.time[-1])
                        hist, _ = nsb_from_dc_histogram(night.dc_report, dc_key)
                        row.update({"dc_scale": scale, "nsb_dc_bkg_dc": bkg, "nsb_dc_hist_dc": hist,
                                    "nsb_dc_bkg": bkg * scale, "nsb_dc_hist": hist * scale})
                    except (ValueError, KeyError, ZeroDivisionError):
                        pass
                if dark:
                    row["nsb_default"], row["nsb_source"] = row["nsb_dark"], "dark"
                elif dc_key == "lst1":
                    row["nsb_default"], row["nsb_source"] = row.get("nsb_dc_hist", np.nan), "dc_hist"
                else:
                    row["nsb_default"], row["nsb_source"] = row.get("nsb_dc_bkg", np.nan), "dc_bkg"
                row["star_mean"] = float(np.mean(good_time_lightcurve(lc, obs.good_time).signal))
                rows.append(row)
    return pd.DataFrame(rows).set_index(["night", "observation", "telescope"])


def measured_snr_table(log, nights=None, n_bins: Sequence[int] = N_BINS, nsb: Optional[pd.DataFrame] = None,
                       lcs: Optional[LightCurves] = None, detrend_window_s: float = DETREND_WINDOW_S,
                       signal_band_hz: Tuple[float, float] = SIGNAL_BAND_HZ, *, include_excluded=False) -> pd.DataFrame:
    """
    Measured SNR of every star observation, telescope and bin size (index
    night, observation, telescope, n_bin), one column per estimator (see the
    module docstring), with ``nsb`` from ``nsb_table``. Also per row (independent of n_bin):
    ``inband_rms`` (fractional-flux noise rms inside ``signal_band_hz``, gap-aware Welch),
    ``mf_sigma`` (noise seen by the reference diffraction template, spectra.matched_filter_sigma),
    ``welch_short`` (1 if a run was too short for a 2**16-sample Welch segment).
    """
    lcs = _lcs(log, nights, lcs, include_excluded)
    # If no nsb table is supplied, we calculate it in the next step
    nsb = nsb if nsb is not None else nsb_table(log, nights, n_bins[0], lcs, include_excluded=include_excluded)
    rows = []
    for night in log.select(nights, include_excluded):
        for ch in night.channels:
            dc_key = night.telescope(ch).dc_key
            tel_name = night.telescope(ch).name
            for obs in night.stars.values():
                if not obs.has_cache:
                    continue
                n = nsb.loc[(night.name, obs.name, tel_name)]
                lc_full = lcs[(night.name, obs.name, ch)]
                lc = good_time_lightcurve(lc_full, obs.good_time)
                comb = despiked = dctrend = None
                bad_frac = np.nan
                if _has_dc(night):
                    try:
                        comb = combine_lightcurve_with_dc(lc, night.dc_report, DC_CUT_HZ[dc_key], telescope=dc_key)
                        despiked, bad_frac = despike_lightcurve_with_dc(lc_full, night.dc_report, telescope=dc_key)
                        despiked = good_time_lightcurve(despiked, obs.good_time)
                        dctrend = good_time_lightcurve(
                            detrend_lightcurve_with_dc(lc_full, night.dc_report, telescope=dc_key), obs.good_time)
                    except (ValueError, KeyError):
                        comb = despiked = dctrend = None
                        bad_frac = np.nan
                nsb_dc = n["nsb_default"] / n["dc_scale"]
                nsb0 = n["nsb_default"]
                nsb_err = n["nsb_dark_err"] if n["nsb_source"] == "dark" else 0.0
                try:
                    band_noise = noise_band_numbers(lc, signal_band_hz, default_template_residual(round(lc.dt, 9)))
                except ValueError:       # no run long enough for a Welch estimate
                    band_noise = {"inband_rms": np.nan, "welch_short": np.nan, "mf_sigma": np.nan}
                for n_bin in n_bins:
                    row = {"night": night.name, "observation": obs.name, "telescope": tel_name,
                           "n_bin": n_bin, "bin_ms": n_bin * lc.dt * 1e3, "nsb_source": n["nsb_source"],
                           "dc_bad_frac": bad_frac, **band_noise,
                           "mean_fast": measured_snr(lc, n_bin, nsb0, nsb_err)["snr"],
                           "peak_fast": measured_snr(lc, n_bin, nsb0, nsb_err, signal="peak")["snr"],
                           "mean_fast_dcnsb": measured_snr(lc, n_bin, n.get("nsb_dc_bkg", np.nan))["snr"],
                           "detr": detrended_snr(lc, n_bin, nsb0, detrend_window_s)["snr"],
                           "detr_full": detrended_snr(lc_full, n_bin, nsb0, detrend_window_s)["snr"]}
                    if comb is not None:
                        row.update({"mean_comb": measured_snr(comb, n_bin, nsb_dc, nsb_err / n["dc_scale"])["snr"],
                                    "peak_comb": measured_snr(comb, n_bin, nsb_dc, nsb_err / n["dc_scale"], signal="peak")["snr"],
                                    "detr_despiked": detrended_snr(despiked, n_bin, nsb0, detrend_window_s)["snr"],
                                    "detr_dcdetrend": detrended_snr(dctrend, n_bin, nsb0, detrend_window_s)["snr"]})
                    rows.append(row)
    return pd.DataFrame(rows).set_index(["night", "observation", "telescope", "n_bin"])


_INSTRUMENTS = {}


def instrument_for(telescope_name: str, extinction_k_B: Optional[float] = None):
    """
    The (cached) Instrument model of a telescope: LST for LST-*, MAGIC otherwise, at La Palma.
    ``extinction_k_B`` (mag/airmass, e.g. the observation log's) rescales the MODTRAN extinction.
    """
    from ..simulation.instruments import Instrument

    kind = "LST" if telescope_name.startswith("LST") else "MAGIC"
    if kind not in _INSTRUMENTS:
        _INSTRUMENTS[kind] = Instrument(kind, site="LaPalma")
    _INSTRUMENTS[kind].set_extinction_k_B(extinction_k_B)
    return _INSTRUMENTS[kind]


def logged_psf_containment(instr, psf_mm) -> float:
    """
    Star light inside the pixel from the PSF logged in elog_notes [mm]: the MAGIC value is the 39 %
    containment radius ('r39', sigma = PSF), the LST-1 one is D80. 1.0 when nothing was logged.
    """
    return instr.psf_containment(psf_mm, "d80" if instr.telescope_type == "LST" else "r39")


def theoretical_snr_table(measured: pd.DataFrame, log, conditions: pd.DataFrame,
                          bands: Sequence[str] = BANDS) -> pd.DataFrame:
    """
    Theoretical SNR for every row of ``measured_snr_table`` (same index):
    theo_shot/zenith/obs_<band> (see the module docstring), plus the airmass
    and mirror_throughput used for theo_obs. The faulty panels are the
    observation's own ``faulty_mirrors`` (the night's value unless overridden).
    """
    import astropy.units as u

    rows = []
    for (night_name, obs_name, tel_name, n_bin), row in measured.iterrows():
        night = log.nights[night_name]
        obs = night.observations[obs_name]
        instr = instrument_for(tel_name, (getattr(log, "site", None) or {}).get("extinction_k", {}).get("B"))
        dc_key = next(night.telescope(ch).dc_key for ch in night.channels
                      if night.telescope(ch).name == tel_name)
        psf = logged_psf_containment(instr, (night.elog_notes or {}).get(dc_key))
        throughput = throughput_from_faulty_panels(obs.faulty_mirrors.get(dc_key), dc_key)
        airmass = conditions["airmass_mean"].get((night_name, obs_name), np.nan) \
            if "airmass_mean" in conditions else np.nan
        tb = row["bin_ms"] * u.ms
        teff = teff_from_spectral_type(obs.spectral_type)
        out = {"night": night_name, "observation": obs_name, "telescope": tel_name, "n_bin": n_bin,
               "airmass": airmass, "mirror_throughput": throughput,
               "teff_K": np.nan if teff is None else teff, "psf_containment": psf}
        for band in bands:
            mag = obs.magnitude(band)
            if mag is None:
                continue
            out[f"theo_shot_{band}"] = instr.signal_to_noise_ratio(mag, tb, temperature_K=teff)
            out[f"theo_zenith_{band}"] = instr.realistic_signal_to_noise_ratio(mag, tb, temperature_K=teff, psf_containment=psf)
            if np.isfinite(airmass):
                out[f"theo_obs_{band}"] = instr.realistic_signal_to_noise_ratio(
                    mag, tb, airmass=airmass, extinction=True, throughput_factor=throughput,
                    temperature_K=teff, psf_containment=psf)
        rows.append(out)
    return pd.DataFrame(rows).set_index(["night", "observation", "telescope", "n_bin"])


def snr_report(log, nights=None, n_bins: Sequence[int] = N_BINS, bands: Sequence[str] = BANDS,
               detrend_window_s: float = DETREND_WINDOW_S, lcs: Optional[LightCurves] = None,
               signal_band_hz: Tuple[float, float] = SIGNAL_BAND_HZ, *,
               include_excluded=False) -> pd.DataFrame:
    """
    Measured and theoretical SNR of every star observation of the selected
    nights (None = all not excluded, see ObservationLog.select), joined on (night, observation, telescope, n_bin), with
    ``ratio_detr_B`` / ``ratio_mean_B`` = detr / mean_fast over theo_obs_B.
    """
    lcs = _lcs(log, nights, lcs, include_excluded)
    conditions = conditions_table(log, nights, include_excluded)
    nsb = nsb_table(log, nights, n_bins[0], lcs, include_excluded=include_excluded)
    measured = measured_snr_table(log, nights, n_bins, nsb, lcs, detrend_window_s, signal_band_hz,
                                  include_excluded=include_excluded)
    theory = theoretical_snr_table(measured, log, conditions, bands)
    results = measured.join(theory)
    if "theo_obs_B" in results:
        results["ratio_detr_B"] = results["detr"] / results["theo_obs_B"]
        results["ratio_mean_B"] = results["mean_fast"] / results["theo_obs_B"]
    return results


def plot_snr_vs_magnitude(results: pd.DataFrame, log, n_bin: int,
                          snr_types: Union[str, Sequence[str]] = ("detr", "mean_fast"),
                          bands: Union[str, Sequence[str]] = ("B",)):
    """
    Measured SNR against magnitude, one panel per band with all telescopes
    together (colour = night, marker = telescope, fill style = SNR type).

    ``snr_types`` are columns of ``results`` (``detr``, ``mean_fast``,
    ``peak_comb``, ...) or theory prefixes (``theo_obs``, ``theo_zenith``,
    ``theo_shot``, resolved to ``<prefix>_<band>``). Returns the figure.
    """
    import astropy.units as u
    import matplotlib.pyplot as plt

    snr_types = (snr_types,) if isinstance(snr_types, str) else tuple(snr_types)
    bands = (bands,) if isinstance(bands, str) else tuple(bands)
    sub = results.xs(n_bin, level="n_bin")

    def column(t, band):
        for c in (f"{t}_{band}", t):
            if c in sub:
                return c
        raise ValueError(f"unknown SNR type {t!r}; available columns: {list(sub.columns)}")

    bin_ms = sub["bin_ms"].iloc[0]
    night_names = list(dict.fromkeys(sub.index.get_level_values("night")))
    colors = dict(zip(night_names, plt.rcParams["axes.prop_cycle"].by_key()["color"]))
    markers = {"MAGIC-1": "s", "MAGIC-2": "o", "LST-1": "^"}
    styles = [{}, {"fill": False}, {"marker": "x"}, {"marker": "+"}]   # per SNR type
    models = {"MAGIC": ("MAGIC-1", "grey"), "LST": ("LST-1", "k")}
    mag_grid = np.linspace(5.5, 9.0, 40)

    def scatter(ax, x, y, tel_name, c, style):
        m = style.get("marker", markers.get(tel_name, "o"))
        if style.get("fill", True):
            ax.scatter(x, y, marker=m, color=c, s=50)
        else:
            ax.scatter(x, y, marker=m, facecolor="none", edgecolor=c, s=50)

    fig, axes = plt.subplots(len(bands), 1, figsize=(9, 5.5 * len(bands)), squeeze=False)
    for ax, band in zip(axes[:, 0], bands):
        cols = [column(t, band) for t in snr_types]
        for tel_example, lc in models.values():
            instr = instrument_for(tel_example)
            ax.plot(mag_grid, [instr.realistic_signal_to_noise_ratio(m, bin_ms * u.ms) for m in mag_grid],
                    "-", color=lc)
        for (night_name, obs_name, tel_name), r in sub.iterrows():
            mag = log.nights[night_name].observations[obs_name].magnitude(band)
            if mag is None:
                continue
            c = colors[night_name]
            for col, style in zip(cols, styles):
                # The actual plotting of the snr that we've selected
                scatter(ax, mag, r[col], tel_name, c, style)
            theo = f"theo_obs_{band}"
            #if show_theory and theo in r and theo not in cols:
                # This is where we plot the most realistic theoretical SNR, with faulty mirrors and all
            ax.scatter(mag, r[theo], marker="_", color=c, s=200)

        ax.set_yscale("log"); ax.invert_xaxis(); ax.grid(alpha=0.3, which="both")
        ax.set_xlabel(f"{band} magnitude"); ax.set_ylabel("SNR")
        ax.set_title(f"{band} band, {bin_ms:.1f} ms bins")
    handles = [plt.Line2D([], [], color=c, marker="o", ls="", label=nm) for nm, c in colors.items()]
    handles += [plt.Line2D([], [], color="grey", marker=m, ls="", label=t) for t, m in markers.items()]
    for t, style in zip(snr_types, styles):
        fill = "grey" if style.get("fill", True) else "none"
        handles.append(plt.Line2D([], [], color="grey", marker=style.get("marker", "o"), mfc=fill, ls="", label=t))
    #if show_theory and "theo_obs" not in snr_types:
    handles.append(plt.Line2D([], [], color="grey", marker="_", ms=15, ls="",
                                label="theory at obs. airmass + mirrors"))
    for name, (_, lc) in models.items():
        handles += [plt.Line2D([], [], color=lc, ls="-", label=f"{name} theory, zenith")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8)
    fig.tight_layout(rect=(0, 0.08 / len(bands) + 0.04, 1, 1))
    return fig
