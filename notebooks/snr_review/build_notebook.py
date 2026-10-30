"""
Build snr_theory_vs_measured.ipynb: a notebook comparing the measured SNR of the test-star
observations with the theoretical one, for review outside the package (it imports only
Instrument, for the star photon counts).

The observation table is taken from the observation log at build time and
embedded in the notebook; everything else is code ported from the package
(see the comments in each cell). Rebuild after editing the log:

    python notebooks/snr_review/build_notebook.py
"""

import os
import pprint

import nbformat as nbf

from kbo_occultation.observations.log import load_observation_log

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "snr_theory_vs_measured.ipynb")


def observation_table(log):
    nights = {}
    for night in log.select():
        obs = {}
        for o in night.observations.values():
            obs[o.name] = {
                "file": os.path.basename(o.path),
                "type": o.type,
                "mag_B": o.mag_B,
                "ra_deg": o.ra_deg,
                "dec_deg": o.dec_deg,
                "good_time": o.good_time,
                "faulty_mirrors": {k: v for k, v in o.faulty_mirrors.items() if v is not None},
            }
        dc = None
        if night.dc_report:
            dc = "/".join(os.path.normpath(night.dc_report).split(os.sep)[-3:])
        nights[night.name] = {"channels": list(night.channels), "dark_run": night.dark_run,
                              "dc_report": dc, "observations": obs}
    return nights


def md(text):
    return nbf.v4.new_markdown_cell(text.strip("\n"))


def code(text):
    return nbf.v4.new_code_cell(text.strip("\n"))


def build():
    log = load_observation_log()
    nights = pprint.pformat(observation_table(log), sort_dicts=False, width=110)
    k_b = log.site["extinction_k"]["B"]

    cells = [
        md("""
# Theoretical vs measured SNR of the fast-photometry test stars

The measured signal-to-noise ratio of the test-star observations is far below the
theoretical one. This notebook shows every input that enters the comparison:

1. the observations and the conditions that enter the theory (magnitude, airmass,
   extinction, faulty mirror panels) and the NSB subtracted from the measurement;
2. the theoretical SNR model, with all corrections (airmass extinction and
   scintillation, mirror area lost to faulty panels);
3. the DC reports of every pixel over each whole night;
4. the light curves (fast signal averaged to 100 ms) with the DC of the on-target pixel;
5. the measured SNR, three ways (mean, detrended, fast+DC combined), against the theory.

All code is included here (ported from the `kbo_occultation` package), so only the
data files are needed.

**Files** (paths relative to `DATA_DIR`, set in the next cell):
```
observations/*.npz             one per run: outlier-cleaned variance per DAQ channel
DCs/<YYYY_MM_DD>/dc_report.pkl slow DC reports (magic_spysii DcReportReader pickles)
reference/                     instrument/atmosphere curves of the theoretical model:
    LST_ref_200_1100_190211a.dat                mirror reflectivity
    LST_qe_R11920-100-02.dat                    PMT quantum efficiency
    MST-N_No7-10_ave.dat                        MAGIC camera window transmission
    LST_transmission_lst_window_No7-10_ave.dat  LST camera window transmission
    atm_trans_2147_1_3_0_0_0.dat                La Palma atmospheric optical depth
    Spectra_NSB_ref.txt                         night-sky background spectrum
```

**Requirements:** numpy, scipy, pandas, matplotlib, astropy, and `magic_spysii`
(only to unpickle the DC reports).

**Flux proxy.** The DAQ records, every 2¹⁸ digitizations of 2 ns (≈ 0.52 ms), the
variance (std²) of the AC-coupled PMT signal, which is proportional to the photon
rate. The star signal is this variance minus the NSB + electronic-pedestal level.
"""),
        code(f"""
import os
import pickle
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.io import ascii
from astropy.time import Time
from scipy import interpolate
from scipy.ndimage import median_filter

from kbo_occultation.simulation.instruments import Instrument  # star photon counts

warnings.filterwarnings("ignore", category=RuntimeWarning)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
pd.set_option("display.precision", 3)

DATA_DIR = "."                     # folder holding observations/, DCs/ and reference/

# DAQ channel -> (telescope, key of the telescope in the DC report)
CHANNELS = {{"A": ("MAGIC-2", "m2"), "B": ("LST-1", "lst1"), "C": ("MAGIC-1", "m1")}}
ORM = EarthLocation(lat=28.7619 * u.deg, lon=-17.8900 * u.deg, height=2200.0 * u.m)

K_EXT_B = {k_b}                     # B-band extinction coefficient at the ORM [mag / airmass]
N_PANELS = {{"m1": 247, "m2": 247, "lst1": 198}}   # mirror panels per telescope
N_BINS = (10, 200)                 # samples per SNR bin: ~5 ms (search cadence) and ~100 ms
DETREND_WINDOW_S = 10.0            # rolling-median baseline of the detrended SNR
DC_CUT_HZ = {{"m1": 0.3, "m2": 0.3, "lst1": 0.15}}  # fast/DC crossover (config.DC_CUT_HZ)

# From the observation log. good_time: [start, stop] s after the first sample (None = whole
# run); faulty_mirrors: {{telescope: panels not correctly focused}} (AMC report, median over the run).
NIGHTS = {nights}
"""),
        md("""
## 1. Loaders and helpers
"""),
        code("""
class LightCurve:
    def __init__(self, time, signal):
        self.time, self.signal = time, signal
        self.dt = float(time[1] - time[0])


def obs_path(name):
    return os.path.join(DATA_DIR, "observations", name)


def load_lightcurve(npz_file, channel):
    \"\"\"Outlier-cleaned variance of one channel; the time grid is rebuilt at the fixed cadence.\"\"\"
    with np.load(obs_path(npz_file), allow_pickle=False) as npz:
        t0, dt, n = float(npz["t0"]), float(npz["dt"]), int(npz["n"])
        signal = np.asarray(npz[f"var_{channel}_clean"], dtype=float)
    return LightCurve(t0 + np.arange(n) * dt, signal)


def good_time_lightcurve(lc, window_rel_s):
    if window_rel_s is None:
        return lc
    t_rel = lc.time - lc.time[0]
    mask = (t_rel >= window_rel_s[0]) & (t_rel <= window_rel_s[1])
    return LightCurve(lc.time[mask], lc.signal[mask])


_DC_REPORTS = {}


def dc_report(night_name):
    \"\"\"The night's DC report (DcReportReader), loaded once; None if the night has none.\"\"\"
    rel = NIGHTS[night_name]["dc_report"]
    if rel is None or not os.path.exists(os.path.join(DATA_DIR, rel)):
        return None
    if night_name not in _DC_REPORTS:
        with open(os.path.join(DATA_DIR, rel), "rb") as f:
            _DC_REPORTS[night_name] = pickle.load(f)
    return _DC_REPORTS[night_name]


def dc_series(night_name, dc_key, pixel=0):
    \"\"\"(unix time, DC) of one pixel: 0 = on target, 1 (and 2 for MAGIC) = background.\"\"\"
    report = dc_report(night_name)
    if report is None:
        raise ValueError(f"no DC report for {night_name}")
    data = getattr(report, dc_key)
    if len(data["times"]) == 0:
        raise ValueError(f"no DC data for {dc_key!r} on {night_name}")
    return Time(data["times"], format="isot", scale="utc").unix, np.asarray(data["dcs"][:, pixel])


def bin_average(t, signal, n):
    \"\"\"Average non-overlapping chunks of n samples; returns (t, mean, std/sqrt(n)) per bin.\"\"\"
    n_trim = len(signal) - len(signal) % n
    tc, sc = t[:n_trim].reshape(-1, n), signal[:n_trim].reshape(-1, n)
    return tc.mean(axis=1), sc.mean(axis=1), sc.std(axis=1) / np.sqrt(n)


def robust_sigma(x):
    x = np.asarray(x)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


lcs = {(night_name, obs_name, ch): load_lightcurve(o["file"], ch)
       for night_name, night in NIGHTS.items() for obs_name, o in night["observations"].items()
       if os.path.exists(obs_path(o["file"])) for ch in night["channels"]}


def stars(night_name):
    return {k: o for k, o in NIGHTS[night_name]["observations"].items() if o["type"] == "star"}
"""),
        md("""
## 2. Observations and the inputs of the SNR comparison

Per star and telescope:
- **Zd / airmass**: mean over the run, computed with astropy from the coordinates and the run times.
- **faulty panels → mirror throughput**: fraction of mirror area left, `1 − n_faulty / n_panels`
  (247 panels for MAGIC, 198 for LST-1; LST-1's AMC is not recorded, so its throughput is 1).
- **extinction factor**: `10^(0.4 k_B (X − 1))`, k_B = 0.23. It is only listed here: the theory
  applies the extinction itself, through the La Palma optical-depth curve at the run's airmass.
- **NSB** (+ electronic pedestal, in fast variance units), subtracted from the measured signal:
  - `dark`: mean of the night's dark run, binned like the stars (2025-12-16, 2026-08-19);
  - `dc_bkg`: the lower of the DC background pixels during the run (MAGIC on 2026-03-08, no dark run);
  - `dc_hist`: first peak of the on-target pixel's DC histogram over the night (LST-1 on 2026-03-08,
    whose background pixel reads too high).

  The DC estimates are converted to fast units with `dc_scale` = mean(fast) / mean(on-target DC) over the run.
"""),
        code("""
def pointing(ra_deg, dec_deg, t_start, t_stop, n_points=20):
    times = Time(np.linspace(t_start, t_stop, n_points), format="unix", scale="utc")
    alt = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg).transform_to(
        AltAz(obstime=times, location=ORM)).alt.deg
    zd = 90.0 - alt
    return float(zd.mean()), float(np.mean(1.0 / np.cos(np.radians(zd))))


def throughput_from_faulty_panels(n_faulty, dc_key):
    return 1.0 if n_faulty is None else 1.0 - n_faulty / N_PANELS[dc_key]


def extinction_factor(airmass, k_mag):
    return 10 ** (0.4 * k_mag * (airmass - 1.0))


def histogram_peak(data, which="last", n_bins=200, min_count=20):
    \"\"\"First local maximum of the histogram, scanning from the high (last) or low (first) end.\"\"\"
    counts, edges = np.histogram(np.asarray(data), bins=n_bins)
    centres = 0.5 * (edges[1:] + edges[:-1])
    counts = np.where(counts < min_count, 0, counts)
    order = range(len(counts) - 1, -1, -1) if which == "last" else range(len(counts))
    best_count, best_idx = -1, None
    for i in order:
        if counts[i] >= best_count:
            best_count, best_idx = counts[i], i
        else:
            return float(centres[best_idx])
    return float("nan")


def nsb_from_dark(lc_dark, n_bin):
    _, signal, sigma = bin_average(lc_dark.time, lc_dark.signal, n_bin)
    return float(np.mean(signal)), float(np.mean(sigma))


def nsb_from_dc_background(night_name, dc_key, t_start, t_stop):
    \"\"\"Lower of the background pixels' mean DC during the run (DC units).\"\"\"
    levels = []
    for pixel in (1, 2):
        try:
            t, v = dc_series(night_name, dc_key, pixel)
        except IndexError:          # LST-1 has a single background pixel
            break
        levels.append(v[(t >= t_start) & (t <= t_stop)])
    return float(np.mean(min(levels, key=np.mean)))


def nsb_from_dc_histogram(night_name, dc_key):
    _, v = dc_series(night_name, dc_key, 0)
    return histogram_peak(v, which="first", n_bins=100, min_count=20)


def dc_to_fast_scale(lc, night_name, dc_key):
    t, v = dc_series(night_name, dc_key, 0)
    sel = (t >= lc.time[0]) & (t <= lc.time[-1])
    return float(np.mean(lc.signal) / np.mean(v[sel]))


rows = []
for night_name, night in NIGHTS.items():
    dark_name = night["dark_run"]
    for ch in night["channels"]:
        tel, dc_key = CHANNELS[ch]
        dark = None
        if dark_name and (night_name, dark_name, ch) in lcs:
            dark = nsb_from_dark(good_time_lightcurve(lcs[(night_name, dark_name, ch)],
                                                      night["observations"][dark_name]["good_time"]), N_BINS[0])
        for obs_name, o in stars(night_name).items():
            if (night_name, obs_name, ch) not in lcs:
                continue
            lc = lcs[(night_name, obs_name, ch)]
            zd, airmass = pointing(o["ra_deg"], o["dec_deg"], lc.time[0], lc.time[-1])
            n_faulty = o["faulty_mirrors"].get(dc_key)
            row = {"night": night_name, "observation": obs_name, "telescope": tel, "dc_key": dc_key,
                   "mag_B": o["mag_B"], "zd_deg": zd, "airmass": airmass, "faulty_panels": n_faulty,
                   "mirror_throughput": throughput_from_faulty_panels(n_faulty, dc_key),
                   "ext_factor_B": extinction_factor(airmass, K_EXT_B), "good_time": o["good_time"],
                   "nsb_dark": dark[0] if dark else np.nan, "nsb_dark_err": dark[1] if dark else 0.0,
                   "dc_scale": np.nan, "nsb_dc_bkg": np.nan, "nsb_dc_hist": np.nan}
            try:
                scale = dc_to_fast_scale(lc, night_name, dc_key)
                row.update({"dc_scale": scale,
                            "nsb_dc_bkg": nsb_from_dc_background(night_name, dc_key, lc.time[0], lc.time[-1]) * scale,
                            "nsb_dc_hist": nsb_from_dc_histogram(night_name, dc_key) * scale})
            except (ValueError, KeyError, ZeroDivisionError):
                pass
            if dark:
                row["nsb"], row["nsb_source"] = row["nsb_dark"], "dark"
            elif dc_key == "lst1":
                row["nsb"], row["nsb_source"] = row["nsb_dc_hist"], "dc_hist"
            else:
                row["nsb"], row["nsb_source"] = row["nsb_dc_bkg"], "dc_bkg"
            row["star_level"] = float(np.mean(good_time_lightcurve(lc, o["good_time"]).signal))
            rows.append(row)
obs_table = pd.DataFrame(rows).set_index(["night", "observation", "telescope"])
obs_table[["mag_B", "zd_deg", "airmass", "faulty_panels", "mirror_throughput", "ext_factor_B", "good_time",
           "nsb", "nsb_source", "star_level", "dc_scale"]]
"""),
        md("""
## 3. Theoretical SNR model

Per-sample fractional noise of the star signal S (photoelectrons per bin), three terms in quadrature:

$$\\sigma^2 = \\underbrace{F^2\\,\\frac{S + B}{S^2}}_{\\text{shot}} + \\underbrace{\\frac{\\sigma_e^2}{S^2}}_{\\text{electronic}} + \\underbrace{10^{-5}\\,D^{-4/3}\\,X^3\\,e^{-2h/H}\\,/\\,(2t)}_{\\text{scintillation (Young 1967)}}, \\qquad \\mathrm{SNR} = 1/\\sigma$$

- **S** = photo-electrons of a star of Johnson B magnitude m: CALSPEC Vega spectrum × response, scaled by 10^(−0.4(m − 0.03)),
  with response T(λ) = atmosphere(X) × mirror reflectivity × QE × camera window × telescope × funnel transmission.
- **B**: NSB photons in one pixel (dark-sky reference spectrum × telescope transmission × pixel solid angle × mirror area).
- **Corrections at the observation's conditions:** the atmosphere is evaluated at the run's airmass X
  (extinction), the scintillation term uses the same X, and S and B are multiplied by the mirror throughput
  left by the faulty panels.
- Parameters (MAGIC / LST): mirror area 240 / 380.3 m²; pixel FoV 0.1° / 0.102°; telescope transmission 0.8928 / 0.9606;
  funnel 0.85 / 0.892; reflectivity and QE files are LST's for both, scaled by 0.82 and 0.86 for MAGIC;
  F = 1.15, σ_e = 0.7 p.e.; D = 2√(area/π); site altitude h = 2200 m, H = 8000 m.
- **Not in the model:** the statistical error of the variance estimate itself (√(2k/2¹⁸) per 0.5 ms frame,
  independent of the star brightness), clouds/transparency, star centring in the pixel, and correlated noise from
  tracking, wind or mirror vibrations.
"""),
        code("""
REF_DIR = os.path.join(DATA_DIR, "reference")
TELESCOPE = {
    "MAGIC": {"area": 240.0, "pixel_fov": 0.1, "transmission": 0.8928, "funnel": 0.85,
              "reflectivity_scale": 0.82, "qe_scale": 0.86, "camera": "MST-N_No7-10_ave.dat", "camera_percent": True},
    "LST": {"area": 380.271, "pixel_fov": 0.102, "transmission": 0.9606, "funnel": 0.892,
            "reflectivity_scale": 1.0, "qe_scale": 1.0, "camera": "LST_transmission_lst_window_No7-10_ave.dat",
            "camera_percent": False},
}
EXCESS_NOISE_F = 1.15
ELEC_NOISE_PE = 0.7
SITE_ALTITUDE_M = 2200.0
WAVELENGTHS = np.arange(100, 1000, 1)


def _curve(filename, wavelength, scale=1.0):
    data = ascii.read(os.path.join(REF_DIR, filename), data_start=1, guess=False, delimiter=",",
                      format="no_header")
    f = interpolate.interp1d(data["col1"], data["col2"])
    out = np.zeros_like(wavelength, dtype=float)
    mask = (wavelength >= np.min(data["col1"])) & (wavelength <= np.max(data["col1"]))
    out[mask] = f(wavelength[mask])
    return out * scale


def telescope_transmission(kind, wavelength):
    p = TELESCOPE[kind]
    return (_curve("LST_ref_200_1100_190211a.dat", wavelength, p["reflectivity_scale"])
            * _curve("LST_qe_R11920-100-02.dat", wavelength, p["qe_scale"])
            * p["transmission"] * p["funnel"]
            * _curve(p["camera"], wavelength, 0.01 if p["camera_percent"] else 1.0))


_ATM_COLS = ["wavelength(nm)"] + [f"c{i}" for i in range(39)] + ["100.000"]


def atmospheric_transmission(wavelength, airmass):
    \"\"\"exp(-X tau(lambda)), tau = optical depth from the site (2.147 km) up to 100 km (last column).\"\"\"
    data = ascii.read(os.path.join(REF_DIR, "atm_trans_2147_1_3_0_0_0.dat"), format="no_header",
                      names=_ATM_COLS, data_start=0, guess=False, delimiter=" ")
    f = interpolate.interp1d(data["wavelength(nm)"], data["100.000"])
    out = np.zeros_like(wavelength, dtype=float)
    mask = (wavelength >= np.min(data["wavelength(nm)"])) & (wavelength <= np.max(data["wavelength(nm)"]))
    out[mask] = np.exp(-airmass * f(wavelength[mask]))
    return out


def star_photons(kind, magnitude, bin_ms, airmass):
    # CALSPEC Vega x response, from the package (La Palma atmosphere = atm_trans_2147 file)
    inst = Instrument(kind, site="LaPalma")
    inst.set_airmass(airmass)
    return inst.num_photons_from_mag(bin_ms, magnitude)


def nsb_photons(kind, bin_ms):
    nsb = ascii.read(os.path.join(REF_DIR, "Spectra_NSB_ref.txt"), format="no_header",
                     names=["wavelength(nm)", "flux"], data_start=0, guess=False, delimiter=" ")
    per_ns = np.sum(nsb["flux"] * telescope_transmission(kind, nsb["wavelength(nm)"]))  # 1/(m2 ns sr)
    pixel_sr = np.pi * (TELESCOPE[kind]["pixel_fov"] / 2) ** 2 * (np.pi / 180) ** 2
    return per_ns * TELESCOPE[kind]["area"] * pixel_sr * bin_ms * 1e6


def noise_budget(kind, magnitude, bin_ms, airmass=1.0, throughput=1.0):
    S = star_photons(kind, magnitude, bin_ms, airmass) * throughput
    B = nsb_photons(kind, bin_ms) * throughput
    D_cm = 2.0 * np.sqrt(TELESCOPE[kind]["area"] / np.pi) * 100.0
    var_shot = EXCESS_NOISE_F ** 2 * (S + B) / S ** 2
    var_elec = ELEC_NOISE_PE ** 2 / S ** 2
    var_scint = 10e-6 * D_cm ** (-4 / 3) * airmass ** 3 * np.exp(-2 * SITE_ALTITUDE_M / 8000.0) / (2 * bin_ms * 1e-3)
    return {"star_photons": S, "nsb_photons": B, "sigma_shot": np.sqrt(var_shot),
            "sigma_elec": np.sqrt(var_elec), "sigma_scint": np.sqrt(var_scint),
            "sigma_total": np.sqrt(var_shot + var_elec + var_scint)}


def theoretical_snr(telescope, magnitude, bin_ms, airmass=1.0, throughput=1.0):
    kind = "LST" if telescope.startswith("LST") else "MAGIC"
    return 1.0 / noise_budget(kind, magnitude, bin_ms, airmass, throughput)["sigma_total"]


# Noise budget of every observation at its conditions, 5 ms bins
dt = next(iter(lcs.values())).dt
budget = pd.DataFrame({idx: noise_budget("LST" if idx[2].startswith("LST") else "MAGIC", r["mag_B"],
                                         N_BINS[0] * dt * 1e3, r["airmass"], r["mirror_throughput"])
                       for idx, r in obs_table.iterrows()}).T
budget.index.names = obs_table.index.names
budget
"""),
        md("""
## 4. DC reports over the whole night

Every pixel of each telescope's slow DC report (pixel 0 on the target, the others background), over the
full night. Grey bands mark the fast-photometry runs.
"""),
        code("""
DC_TELESCOPES = {"m1": "MAGIC-1", "m2": "MAGIC-2", "lst1": "LST-1"}
to_dt = lambda t: Time(t, format="unix").datetime


def run_spans(night_name):
    return {name: (lcs[(night_name, name, ch)].time[0], lcs[(night_name, name, ch)].time[-1])
            for name in NIGHTS[night_name]["observations"]
            for ch in NIGHTS[night_name]["channels"][:1] if (night_name, name, ch) in lcs}


for night_name in NIGHTS:
    report = dc_report(night_name)
    if report is None:
        continue
    tels = {k: getattr(report, k) for k in DC_TELESCOPES
            if getattr(report, k, None) is not None and len(getattr(report, k)["times"]) > 0}
    spans = run_spans(night_name)
    fig, axes = plt.subplots(len(tels), 1, squeeze=False, figsize=(14, 3 * len(tels)), sharex=True)
    for ax, (key, data) in zip(axes[:, 0], tels.items()):
        t = Time(data["times"], format="isot", scale="utc").datetime
        dcs = np.asarray(data["dcs"])
        for p in range(dcs.shape[1]):
            ax.plot(t, dcs[:, p], lw=0.7, label=f"pixel {p} ({'on-target' if p == 0 else 'background'})")
        for t0, t1 in spans.values():
            ax.axvspan(*to_dt([t0, t1]), color="grey", alpha=0.2)
        ax.set_ylabel(f"{DC_TELESCOPES[key]} DC")
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(alpha=0.3)
    for name, (t0, t1) in spans.items():
        axes[0, 0].text(to_dt(t0), axes[0, 0].get_ylim()[1], name, rotation=90, va="top", fontsize=8)
    axes[-1, 0].set_xlabel("UTC")
    fig.suptitle(f"{night_name}: DC report, all pixels, whole night (grey: fast-photometry runs)")
    fig.tight_layout()
    plt.show()
"""),
        md("""
## 5. Light curves

Fast signal (variance) averaged to 100 ms (blue, left axis) and the DC of the on-target pixel over the same run
(orange, right axis). Stars go from brightest to faintest (B), dark runs last. The green band is the good-time
window used for the SNR (runs without one use the whole run).
"""),
        code("""
def plot_order(night_name):
    st = sorted(stars(night_name).items(), key=lambda kv: (kv[1]["mag_B"] is None, kv[1]["mag_B"] or 0))
    return st + [(k, o) for k, o in NIGHTS[night_name]["observations"].items() if o["type"] != "star"]


for night_name, night in NIGHTS.items():
    order = [(k, o) for k, o in plot_order(night_name) if (night_name, k, night["channels"][0]) in lcs]
    fig, axes = plt.subplots(len(order), len(night["channels"]), squeeze=False,
                             figsize=(5.5 * len(night["channels"]), 2.2 * len(order)), sharex=True)
    for i, (obs_name, o) in enumerate(order):
        label = f"B={o['mag_B']}" if o["type"] == "star" else o["type"]
        for j, ch in enumerate(night["channels"]):
            tel, dc_key = CHANNELS[ch]
            lc = lcs[(night_name, obs_name, ch)]
            t, s, _ = bin_average(lc.time - lc.time[0], lc.signal, int(round(0.1 / lc.dt)))
            ax = axes[i, j]
            ax.plot(t, s, lw=0.8, color="tab:blue")
            ax.set_ylabel("variance", color="tab:blue", fontsize=8)
            ax.tick_params(axis="y", labelcolor="tab:blue")
            try:
                t_dc, v_dc = dc_series(night_name, dc_key, 0)
                sel = (t_dc >= lc.time[0]) & (t_dc <= lc.time[-1])
                if sel.any():
                    ax_dc = ax.twinx()
                    ax_dc.plot(t_dc[sel] - lc.time[0], v_dc[sel], lw=0.8, color="tab:orange")
                    ax_dc.set_ylabel("DC", color="tab:orange", fontsize=8)
                    ax_dc.tick_params(axis="y", labelcolor="tab:orange")
            except ValueError:
                pass
            if o["good_time"]:
                ax.axvspan(*o["good_time"], color="green", alpha=0.2)
            ax.set_title(f"{obs_name} ({label}) – {tel}", fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("time since start [s]")
    fig.suptitle(f"{night_name}: variance (100 ms bins, blue) and DC on-target pixel (orange)")
    fig.tight_layout()
    plt.show()
"""),
        md("""
## 6. Measured SNR

All on the good-time window, at the binned cadence (n_bin samples per point), NSB from section 2 subtracted.

| column | signal | noise |
|---|---|---|
| `mean_fast` | mean of the binned fast signal − NSB | mean within-bin standard error ⊕ NSB error |
| `detr` | rolling-median baseline (10 s) − NSB | robust (MAD) bin-to-bin scatter of signal/baseline − 1, × baseline |
| `mean_comb` | as `mean_fast`, on the fast+DC combined signal (DC units, NSB / dc_scale) | as `mean_fast` |

The combined signal takes the DC report's spectrum below 1 Hz (0.5 Hz for LST-1) and the fast signal's above.
`mean_*` only see the scatter *within* a bin, so noise correlated over one bin or more is invisible to them
(and they run ~8 % high at n_bin = 10 from the small-sample bias of the std). `detr` sees the bin-to-bin
scatter, which is what an occultation dip of that duration has to beat.
"""),
        code("""
def measured_snr(lc, n_bin, nsb, nsb_err=0.0):
    _, binned, sigma = bin_average(lc.time, lc.signal, n_bin)
    star = float(np.mean(binned)) - nsb
    return star / float(np.sqrt(np.mean(sigma) ** 2 + nsb_err ** 2))


def rolling_median_baseline(time, signal, window_s):
    \"\"\"Rolling median over window_s, on block means of window/64 samples, interpolated back.\"\"\"
    dt = float(np.median(np.diff(time)))
    window = max(int(round(window_s / dt)), int(round(window_s / 20 / dt)))
    stride = max(1, window // 64)
    if stride == 1:
        return median_filter(signal, size=window, mode="reflect")
    n_blocks = len(signal) // stride
    sig_blocks = signal[:n_blocks * stride].reshape(-1, stride).mean(axis=1)
    t_blocks = time[:n_blocks * stride].reshape(-1, stride).mean(axis=1)
    return np.interp(time, t_blocks, median_filter(sig_blocks, size=max(window // stride, 3), mode="reflect"))


def detrended_snr(lc, n_bin, nsb, window_s=DETREND_WINDOW_S):
    t, binned, _ = bin_average(lc.time, lc.signal, n_bin)
    baseline = rolling_median_baseline(t, binned, window_s)
    level = float(np.median(baseline))
    return (level - nsb) / (robust_sigma(binned / baseline - 1.0) * level)


def combine_with_dc(lc, night_name, dc_key, cut_hz):
    \"\"\"Fast signal with its spectrum below cut_hz replaced by the DC's (in DC units).\"\"\"
    t_dc, v_dc = dc_series(night_name, dc_key, 0)
    dc_dt = float(np.median(np.diff(t_dc)))
    if lc.time.min() < t_dc.min() - dc_dt or lc.time.max() > t_dc.max() + dc_dt:
        raise ValueError("fast light curve not covered by the DC report")
    dc_on_fast = np.interp(lc.time, t_dc, v_dc)
    freq = np.fft.rfftfreq(len(lc.time), d=float(np.median(np.diff(lc.time))))
    mean_dc, mean_fast = np.mean(dc_on_fast), np.mean(lc.signal)
    spec = np.fft.rfft(dc_on_fast / mean_dc)
    above = freq >= cut_hz
    spec[above] = np.fft.rfft(lc.signal / mean_fast)[above]
    return LightCurve(lc.time, np.fft.irfft(spec, n=len(lc.time)) * mean_dc)


rows = []
for (night_name, obs_name, tel), r in obs_table.iterrows():
    ch = next(c for c in NIGHTS[night_name]["channels"] if CHANNELS[c][0] == tel)
    o = NIGHTS[night_name]["observations"][obs_name]
    lc = good_time_lightcurve(lcs[(night_name, obs_name, ch)], o["good_time"])
    nsb_err = r["nsb_dark_err"] if r["nsb_source"] == "dark" else 0.0
    try:
        comb = combine_with_dc(lc, night_name, r["dc_key"], DC_CUT_HZ[r["dc_key"]])
    except (ValueError, KeyError):
        comb = None
    for n_bin in N_BINS:
        bin_ms = n_bin * lc.dt * 1e3
        rows.append({"night": night_name, "observation": obs_name, "telescope": tel, "n_bin": n_bin,
                     "bin_ms": bin_ms, "mag_B": r["mag_B"],
                     "mean_fast": measured_snr(lc, n_bin, r["nsb"], nsb_err),
                     "detr": detrended_snr(lc, n_bin, r["nsb"]),
                     "mean_comb": measured_snr(comb, n_bin, r["nsb"] / r["dc_scale"]) if comb is not None else np.nan,
                     "theo_obs_B": theoretical_snr(tel, r["mag_B"], bin_ms, r["airmass"], r["mirror_throughput"])})
results = pd.DataFrame(rows).set_index(["night", "observation", "telescope", "n_bin"])
for col in ["mean_fast", "detr", "mean_comb"]:
    results[f"{col}/theo"] = results[col] / results["theo_obs_B"]
results
"""),
        md("""
## 7. Measured vs theoretical SNR

Colour = night, marker = telescope. Filled: detrended; open: mean; ×: fast+DC combined (mean); horizontal tick:
theory at the observation's airmass and mirror throughput. The black line is the same model at the zenith with
all mirrors working, for reference.
"""),
        code("""
night_colors = dict(zip(NIGHTS, plt.rcParams["axes.prop_cycle"].by_key()["color"]))
markers = {"MAGIC-1": "s", "MAGIC-2": "o", "LST-1": "^"}
mag_grid = np.linspace(5.5, 9.0, 40)

for n_bin in N_BINS:
    sub = results.xs(n_bin, level="n_bin")
    bin_ms = sub["bin_ms"].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    for ax, (kind, tel_example) in zip(axes, [("MAGIC", "MAGIC-1"), ("LST", "LST-1")]):
        ax.plot(mag_grid, [theoretical_snr(tel_example, m, bin_ms) for m in mag_grid], "k-", lw=1)
        for (night_name, obs_name, tel), r in sub.iterrows():
            if not tel.startswith(kind):
                continue
            c, m = night_colors[night_name], markers[tel]
            ax.scatter(r["mag_B"], r["detr"], marker=m, color=c, s=50)
            ax.scatter(r["mag_B"], r["mean_fast"], marker=m, facecolor="none", edgecolor=c, s=50)
            ax.scatter(r["mag_B"], r["mean_comb"], marker="x", color=c, s=40)
            ax.scatter(r["mag_B"], r["theo_obs_B"], marker="_", color=c, s=250)
        ax.set_yscale("log"); ax.invert_xaxis(); ax.grid(alpha=0.3, which="both")
        ax.set_xlabel("B magnitude"); ax.set_ylabel(f"SNR in {bin_ms:.1f} ms")
        ax.set_title(f"{kind} telescopes")
    handles = [plt.Line2D([], [], color=c, marker="o", ls="", label=nm) for nm, c in night_colors.items()]
    handles += [plt.Line2D([], [], color="grey", marker=m, ls="", label=t) for t, m in markers.items()]
    handles += [plt.Line2D([], [], color="grey", marker="o", ls="", label="detrended"),
                plt.Line2D([], [], color="grey", marker="o", mfc="none", ls="", label="mean (within-bin)"),
                plt.Line2D([], [], color="grey", marker="x", ls="", label="fast+DC combined"),
                plt.Line2D([], [], color="grey", marker="_", ms=15, ls="", label="theory, obs. conditions"),
                plt.Line2D([], [], color="k", ls="-", label="theory, zenith, all mirrors")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8)
    fig.suptitle(f"Measured vs theoretical SNR, {bin_ms:.1f} ms bins")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    plt.show()

results.groupby(level=["night", "telescope", "n_bin"])[["mean_fast/theo", "detr/theo", "mean_comb/theo"]].median()
"""),
        md("""
## 8. What is known so far

- **Correlated (red) noise dominates.** After a 10 s detrend the bin-to-bin noise is nearly flat at ~1 % from 1 ms
  to 1 s; the power spectrum puts 1–4 % RMS in 0.1–1 Hz against ≲ 0.5 % in 1–10 Hz. Averaging stops helping after
  a few ms, so the 100 ms ratios are much worse than the 5 ms ones. Candidates: wind or mirror vibrations moving
  the star image across the 0.1° pixel, tracking jitter, scintillation above the Young model.
- **Per-frame floor (MAGIC).** Even at 0.5 ms the scatter has a level-independent floor of ≈ 0.4 % (≈ 1.4–1.5 ×
  √(2/2¹⁸)), consistent with the statistical error of the per-frame variance estimate with ~2 correlated
  digitizer samples. It is not in the model. LST-1 shows no such floor.
- **Size of the gap** (medians per night and telescope, table above). At 5 ms the mean estimator reaches
  ~0.1–0.2 of the theory and the detrended one ~0.02–0.08; at 100 ms the detrended one drops to ~0.01–0.025,
  while the mean and combined ones lose less because they don't see the correlated noise. The detrended
  estimator is the relevant one for an occultation search. LST-1 on 2026-08-19 is lowest (mean ≈ 0.03 of theory
  at 5 ms).
- **NSB.** 2026-03-08 has no dark run: its NSB comes from the DC, converted with the fast/DC ratio, which ignores
  the electronic pedestal of the variance.
- **Good-time windows** are not set for 2026-08-19 (whole runs used).
- **Mirrors.** Only the number of panels in error/undefined state enters (as lost area); their actual
  reflectivity/PSF state is not known. LST-1's panel state is not recorded.
"""),
    ]

    nb = nbf.v4.new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"}})
    nbf.write(nb, OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    build()
