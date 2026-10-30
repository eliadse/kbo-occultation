# kbo-occultation

Simulation and detection of stellar occultation diffraction patterns by sub-km
Kuiper Belt Objects (KBOs), for fast photometry with the **MAGIC + LST-1** array at
the Roque de los Muchachos Observatory (ORM).

**Authors:** E. do Souto Espiñeira, T. Hassan, V. Pascual

When a small KBO passes in front of a background star, it casts a Fresnel
diffraction shadow that sweeps across the Earth in milliseconds. The package covers
four tasks:

| I want to… | Start here |
|---|---|
| Add a new observation night | [`kbo-night prepare`](#adding-a-new-night) |
| Measure the SNR of real observations and compare it with theory | [`kbo-night snr`](#measuring-the-snr-of-real-observations) |
| Simulate an occultation light curve or a theoretical SNR | [`compute_lightcurve`, `Instrument`](#simulating-light-curves-and-theoretical-snr) |
| Search a recording for occultations | [`search_observation`](#searching-an-observation-for-occultations) |
| Measure detection efficiency and set upper limits | [injection Monte Carlo](#injection--recovery-and-upper-limits) |

---

## Installation

Requires Python 3.9+ and the scientific-Python stack (numpy, scipy, astropy,
matplotlib, pandas, pyyaml, ruamel.yaml, astroquery — installed automatically).

```bash
git clone <repo-url>
cd kbo_occultation
pip install -e ".[test]"
```

This installs the `kbo_occultation` Python package and the **`kbo-night`** command,
which runs the per-night workflow (`kbo-night --help`). Reading the slow DC reports
(`dc_report.pkl`) additionally needs `magic_spysii` installed, because the pickles
contain its classes.

---

## Package map

```
kbo_occultation/
  config.py           dataclasses (KBO, star, band, grid, numerics) + the telescope array ORM_ARRAY
  cli.py              the kbo-night command (also installed as `kbo`)
  simulation/         forward model: what an occultation looks like, how bright a star is
  observations/       real data: raw files + night metadata -> LightCurves and conditions
  analysis/           measurements on LightCurves: SNR, search, injection, upper limits
  data/               bundled filter / QE / atmosphere / NSB curves, and data/observations/ (your data)
examples/
  simulation/         light curves, sweeps, diffraction PSDs (no data needed)
  real_data/          search and injection Monte Carlo on a logged observation
  studies/2025_12_16/ archived one-off noise studies of the first night
notebooks/
  snr_study.ipynb     the SNR study: kbo-night snr's tables plus noise diagnostics
```

| Module | What it is for |
|---|---|
| `simulation.engine` | `compute_lightcurve`, `OccultationEngine`: polychromatic Fresnel diffraction with finite stellar disk |
| `simulation.sweep` | `run_parameter_sweep` over KBO radius / distance / impact parameter / star size |
| `simulation.instruments` | `Instrument`: telescope throughput, photon counts, `noise_budget`, theoretical SNR |
| `simulation.detectability` | spatial → time conversion, binning, SNR of an aligned dip |
| `simulation.kinematics` | shadow velocity of a KBO (`ShadowKinematics`) |
| `observations.log` | `load_observation_log`: the observation log (nights, stars, good times, mirrors) and each observation's light curves |
| `observations.prepare` | `prepare_night`: writes a night's log entry (filenames, data, SIMBAD, CC/AMC reports) and builds its caches and conditions |
| `observations.simbad` | star coordinates, magnitudes and spectral types from SIMBAD |
| `observations.photometry` | `LightCurve`, `.npz` pre-processing |
| `observations.io` | stat-binary file reader |
| `observations.filtering` | high-pass filters, outlier removal |
| `observations.dc_combine` | slow DC reports: loading, combining with / detrending / despiking the fast signal |
| `observations.cc_reports` | MAGIC CC `.rep` parsing, pointing / moon geometry, `conditions_table` |
| `analysis.snr` | `snr_report`: measured vs theoretical SNR per observation (what `kbo-night snr` runs) |
| `analysis.snr_measurement` | the individual SNR, NSB and noise estimators |
| `analysis.matched_filter`, `analysis.search`, `analysis.coincidence` | blind matched-filter search with χ² veto and multi-telescope coincidence |
| `analysis.bands` | `signal_band`: frequency band of the diffraction signal (f05 … f99, matched-filter band, chord time) |
| `analysis.spectra`, `analysis.surrogates` | gap-aware Welch PSD, in-band noise, noise seen by a matched filter; Gaussian phase-randomised surrogates |
| `analysis.mimic` | candidate catalogue, surrogates, false-alarm rates (what `kbo-night mimic` runs) |
| `analysis.noise_removal` | efficiency at fixed real-data FAR of each noise-removal method (Butterworth scan, PSD whitening, DC detrend / swap / despike veto) |
| `analysis.injection`, `analysis.injection_batch` | signal injection / recovery Monte Carlo |
| `analysis.upper_limits`, `analysis.rates` | detection efficiency → KBO surface-density limits, event rates |

The most used names are also importable from the top level
(`from kbo_occultation import compute_lightcurve, load_observation_log, snr_report, …`).

---

## Data layout

The observation data is not committed. Only the observation log is versioned.

```
kbo_occultation/data/observations/          ($KBO_OBSERVATIONS_DIR overrides this location)
  observation_log.yaml                      the log: one entry per night (versioned)
  Spectrum_stats_..._<name>_..._<UTC>.bin   raw stat-binary recordings (~140 MB each)
  Spectrum_stats_..._<name>_..._<UTC>.npz   compact caches built by `kbo-night prepare`
  DCs/<YYYY_MM_DD>/dc_report.pkl            slow DC reports (optional)
  conditions/<YYYY-MM-DD>/<tel>_<report>.csv  CC-report extracts built by `kbo-night prepare`
<cc_reports_root>/...                       raw MAGIC CC .rep files (~300 MB per telescope and night)
```

`cc_reports_root` is set in the log's `paths:` section and can be overridden with
`$KBO_CC_REPORTS_ROOT`, so the log works on any machine.

About the raw data: the three DAQ channels map to telescopes as **C = MAGIC-1,
A = MAGIC-2, B = LST-1** (`config.ORM_ARRAY`). The per-sample flux proxy is the
**variance** (`std²`) of the 2¹⁸ digitised samples in each record (~0.52 ms), not a
mean count. The `.npz` cache keeps the raw and the outlier-cleaned variance of every
channel plus `t0`/`dt`/`n`.

---

## Adding a new night

1. **Copy the files.** Put the `Spectrum*.bin` files in `kbo_occultation/data/observations/`
   and the DC report in `DCs/<YYYY_MM_DD>/dc_report.pkl`. Put the CC `.rep` files
   (`CC_M1_<YYYY_MM_DD>_*.rep`, `CC_M2_…`) anywhere under `cc_reports_root`.

2. **Prepare the night.**

   ```bash
   kbo-night prepare 2026-10-20
   ```

   This writes the night into `observation_log.yaml` itself; nothing has to be copied
   or pasted. It picks up every `Spectrum*.bin` whose filename timestamp falls between
   12:00 UTC on that date and 12:00 UTC the next day, then:
   - builds the `.npz` caches;
   - condenses the CC reports into `conditions/2026-10-20/`;
   - fills in every field it can derive (table below);
   - prints what it filled and what is still missing.

   Each filled value carries a comment naming its source (`# SIMBAD`, `# CC target`,
   `# AMC`, `# data`, `# found`). It needs network access for SIMBAD; if SIMBAD can't
   be reached, those fields stay `null` and the next run fills them.

3. **Add the manual information (optional).** Edit the night in the log: `good_time`
   windows (`null` = use the whole run), `elog_notes`, `notes`, `optical_filter`
   if a filter was mounted, and `exclude` to leave a problematic night or run out of
   the analysis ([below](#excluding-a-night-or-observation)).

4. **Measure.** `kbo-night snr 2026-10-20` gives the SNR table (next section).
   `kbo-night list` shows every night's status.

| Field | Filled from |
|---|---|
| observations, `file`, `type`, `dark_run`, `daq` | the `.bin` filenames (`type: dark` for names containing "dark") |
| `dc_report` | `DCs/<YYYY_MM_DD>/dc_report.pkl`, if present |
| `cc_reports` | `CC_M<n>_<YYYY_MM_DD>_*.rep` under `cc_reports_root` |
| `channels` | the data: a channel at the electronic pedestal (≈0.4–0.7 variance units, the same in every run) recorded no light and is left out |
| star `ra_deg`, `dec_deg`, `mag_B`, `mag_V`, `spectral_type`, `simbad_id` | SIMBAD, by the name in the filename (`HD-4335` → `HD 4335`); a warning is printed if the position is > 0.5° from the CC target |
| dark-run `ra_deg`, `dec_deg` | the CC-REPORT target RA/Dec during the run |
| `faulty_mirrors` | AMC panels not correctly focused (undefined + error): the night's median, plus a per-run `faulty_mirrors` where a run differs. LST-1 is not in the MAGIC CC reports. |
| `good_time`, `elog_notes`, `notes`, `optical_filter`, `exclude`, `manual`, `manual_conditions` | you (optional) |

Zd, airmass and moon are not stored in the log; they are computed from the
coordinates and the run times.

**Re-running `prepare`** is always safe:
- It only fills empty (`null`) fields. A value you set or corrected by hand is kept.
- Runs that appear later (new `.bin` files of that night) are added to the entry.
- `--refresh` re-derives and overwrites the automatic fields; it never touches the
  manual ones.
- `--dry-run` prints the changes without writing the log.
- `--force` rebuilds the caches and conditions, e.g. after changing the
  outlier-removal parameters.
- `kbo-night prepare` with no date re-checks every night in the log.

### Excluding a night or observation

A night or a single observation that would bias the results can be left out by
giving it an `exclude:` reason in the log (`exclude: true` also works):

```yaml
  "2026-03-08":
    date: "2026-03-08"
    exclude: "no dark run: NSB only from the DC background pixel / DC histogram; biases the SNR studies"
    ...
    observations:
      HD-12345:
        exclude: "tracking lost"     # one run only
```

- `kbo-night snr`, `kbo-night conditions` and the notebooks skip excluded nights and
  observations, and print what they skipped and why.
- A night named explicitly (`kbo-night snr 2026-03-08`) is still used; its excluded
  observations are not.
- `--include-excluded` (CLI) or `include_excluded=True` (Python: `log.select`,
  `snr_report`, `conditions_table`, …) uses everything.
- `kbo-night list` shows the reason; `prepare` and `missing` always work on every
  night.

---

## Measuring the SNR of real observations

```bash
kbo-night snr 2026-10-20 --plot          # one night
kbo-night snr --plot                     # every night in the log, except excluded ones
kbo-night snr --plot --include-excluded  # really every night
kbo-night snr 2026-10-20 --bins 10 50 200 --detrend-window 10
```

These commands write `outputs/snr/snr_<nights>.csv` with all columns, print a summary,
and with `--plot` save SNR-vs-magnitude figures per bin size. The same from Python:

```python
from kbo_occultation import load_observation_log, snr_report

log = load_observation_log()
results = snr_report(log, ["2026-10-20"], n_bins=(10, 200))
results.loc[("2026-10-20", "HD-12345", "MAGIC-1", 10), ["detr", "theo_obs_B", "ratio_detr_B"]]
results = snr_report(log)                          # all nights, except excluded ones
results = snr_report(log, include_excluded=True)   # all nights
```

One row per (night, star observation, telescope, `n_bin`). `n_bin` is the number of
0.52 ms samples per bin (10 ≈ 5 ms, the search cadence). Every SNR uses the good-time
window, with the NSB subtracted.

**Measured SNR**

| Column | Signal | Noise |
|---|---|---|
| **`detr`** | rolling-median baseline − NSB | robust bin-to-bin σ after removing a 10 s rolling baseline |
| `detr_full` | as `detr`, whole run instead of the good-time window | |
| `detr_despiked` | as `detr`, dropping samples near DC excursions (wind / tracking dips) | |
| `detr_dcdetrend` | as `detr`, after dividing by the smoothed DC trend | |
| `mean_fast` / `peak_fast` | mean / last histogram peak of the binned fast signal − NSB | mean within-bin standard error ⊕ NSB error |
| `mean_comb` / `peak_comb` | same on the fast + DC combined signal | |
| `mean_fast_dcnsb` | `mean_fast` with the NSB from the DC background pixel | |

**Use `detr` to compare with theory.** The `mean_*` / `peak_*` estimators only see
scatter *within* a bin. Noise correlated over one bin or longer is invisible to them,
so they run high, by about 8 % more at 5 ms from the small-sample bias.

**NSB** (`nsb_source`): the default is the night's dark run if there is one.
Otherwise it is the first peak of the DC histogram for LST-1 (`dc_hist`) and the
lower DC background pixel for MAGIC (`dc_bkg`), both converted to fast units with
mean(fast) / mean(DC). `dc_hist` is wrong if the camera was closed during part of
the DC report, because the first peak is then the closed-camera level. Take a dark
patch on new nights.

**Theory**, per band `B` / `V`, from `Instrument` (`simulation.instruments`):
- `theo_shot_<band>`: shot noise only;
- `theo_zenith_<band>`: shot × F² + electronic + scintillation, at the zenith;
- `theo_obs_<band>`: the same at the observation's airmass, with extinction, times
  the `mirror_throughput` left after the faulty panels.

`ratio_detr_B = detr / theo_obs_B`.

The individual steps (`conditions_table`, `nsb_table`, `measured_snr_table`,
`theoretical_snr_table`) and the low-level estimators (`analysis.snr_measurement`:
`detrended_snr`, `measured_snr`, `noise_vs_bin_size`, `band_noise_rms`,
`frame_noise`, …) can be called directly.
[notebooks/snr_study.ipynb](notebooks/snr_study.ipynb) uses them for the
conditions plots, noise vs bin size, PSD bands and the per-frame noise floor.
`kbo-night conditions NIGHT [--csv FILE] [--include-excluded]` prints the per-observation conditions on
their own: pointing, moon, weather, pyrometer, starguider and AMC.

Working with a single observation:

```python
obs = log.get("2026-08-19", "HD176914")
lcs = obs.lightcurves(good_time=True)        # {channel: LightCurve}, outlier-cleaned
lc = obs.lightcurve("B", cleaned=False)      # one channel, raw variance
night = log.nights["2026-08-19"]
night.telescope("B").name, night.dark, night.dc_report
```

---

## Simulating light curves and theoretical SNR

No data needed. A diffraction light curve:

```python
import matplotlib.pyplot as plt
from kbo_occultation import (compute_lightcurve, KBOConfig, StarConfig,
                             BandpassConfig, GridConfig, NumericalConfig)

kbo   = KBOConfig(radius_m=500.0, distance_au=40.0)   # a 500 m KBO at 40 AU
star  = StarConfig(temperature_K=20000, angular_radius_mas=0.03)
band  = BandpassConfig(lam_min_nm=400, lam_max_nm=430, n_lambda=25)
grid  = GridConfig(x_max_m=5000, n_x=500)             # shadow-plane sampling
num   = NumericalConfig()                             # numerical resolution

x, intensity = compute_lightcurve(kbo, star, band, grid, num)   # x in metres

plt.plot(x, intensity)
plt.xlabel("Position (m)"); plt.ylabel("Normalized intensity")
plt.show()
```

To scan many KBO radii, distances, impact parameters and star sizes efficiently,
use `run_parameter_sweep` (one `OccultationEngine` per band, Fresnel profiles
cached). `to_dataframe` flattens the result. To get a light curve in time, convert
the spatial axis with `spatial_to_time(x, shadow_velocity_mps)`, where the velocity
comes from `simulation.kinematics.shadow_velocity_opposition` or `ShadowKinematics`.

Theoretical photon counts, noise and SNR of a telescope:

```python
import astropy.units as u
from kbo_occultation import BandpassConfig
from kbo_occultation.simulation.instruments import Instrument

magic = Instrument("MAGIC", site="LaPalma")    # or "LST"; filter_type="MAGIC_SII" (data/optical_filters/)
magic.signal_to_noise_ratio(7.0, 5 * u.ms)     # shot noise only
magic.realistic_signal_to_noise_ratio(7.0, 5 * u.ms, airmass=1.2, extinction=True)
magic.noise_budget(7.0, 5 * u.ms, airmass=1.2) # dict: sigma_shot, sigma_elec, sigma_scint, sigma_floor, sigma_red, sigma_total
band = BandpassConfig(*magic.response_band())  # band covering the instrument response, for compute_lightcurve
```

Examples in [examples/simulation/](examples/simulation/):
- [basic_simulation.py](examples/simulation/basic_simulation.py): monochromatic vs polychromatic vs finite star
- [plot_example.py](examples/simulation/plot_example.py): minimal single curve
- [parameter_sweep_example.py](examples/simulation/parameter_sweep_example.py): `run_parameter_sweep` and a flat-noise detectability cut
- [plot_lightcurve_comparison.py](examples/simulation/plot_lightcurve_comparison.py): light curves at the real cadence with instrument noise
- [plot_diffraction_psd.py](examples/simulation/plot_diffraction_psd.py): power spectrum of the diffraction pattern
- [signal_band_study.py](examples/simulation/signal_band_study.py): convergence and physical scan of the signal band (see [Signal band](#signal-band))

---

### Signal band

Where in frequency does the diffraction signal live? `analysis.bands.signal_band`
answers with several numbers instead of a single cut-off: `f05 … f99` (the
frequency below which that fraction of the power of 1 − I lies, at the shadow
velocity), `f_low = f05`, the matched-filter band `[f05, f95]` and the chord time
`H/v`. Run [signal_band_study.py](examples/simulation/signal_band_study.py) to
reproduce the convergence tests and the scan (`outputs/simulation/signal_band/`).

Reference case: 100 m KBO at 40 AU, 25 km/s, b = 0, 320–450 nm, T = 10 000 K.

| r★ (mas) | f05 | f50 | f90 | f95 | f99 (Hz) | chord time |
|---|---|---|---|---|---|---|
| 0.010 | 0.71 | 13.3 | 22.5 | 26.7 | 31.6 | 0.17 s |
| 0.020 | 0.30 | 4.0 | 14.9 | 16.0 | 21.3 | 0.20 s |
| 0.030 | 0.18 | 2.7 | 6.4 | 8.4 | 21.5 | 0.22 s |
| 0.040 | 0.17 | 2.6 | 6.1 | 7.4 | 15.6 | 0.24 s |

- **Numerics are converged**: f90 and f95 change by < 0.6 % between the two finest
  settings of x_max (30 → 60 km), dx (5 → 2.5 m), n_λ (41 → 81), `n_star_side`
  (128 → 256) and `n_r_grid` (6000 → 12 000). Coarser settings are not safe at large
  r★: `n_star_side = 32` over-estimates f99 of a 0.03 mas star by 9 %.
- **The spectrum is bimodal**: an envelope hump below ~8 Hz (the shadow's width) and
  a Fresnel-fringe hump at 12–16 Hz that carries 22 % of the power at r★ = 0.02 mas,
  10 % at 0.025 mas and 2 % at 0.03 mas. Quantiles are threshold crossings, so
  f90 / f95 hop from the fringe hump to the envelope hump (14.9 → 12.9 → 6.4 Hz for
  r★ = 0.020 → 0.025 → 0.030 mas) once the fringes carry < 10 % / 5 % of the power.
  This also explains the earlier "0.02 mas lower than 0.03 mas" in f99: the 99 %
  point sits in the far tail, which is nearly independent of r★ above 0.02 mas
  (≈ 21 Hz), and was perturbed by coarse disk sampling and by subtracting the mean
  instead of the baseline. Use f90 / f95 for a band, not f99.
- **Scalings**: f ∝ v exactly (f95 = 1.7, 5.1, 8.4, 10.1 Hz at 5, 15, 25, 30 km/s);
  moving from 30 to 50 AU lowers f95 from 16 to 6 Hz; R barely matters for
  R ≲ 200 m (diffraction-dominated; f90 identical for 20–200 m); an impact
  parameter of 1–2 km raises f90 from 6.4 to 8–14 Hz as the chord leaves the geometric
  shadow. The instrument response is a small effect (MAGIC and LST give identical
  bands, they share QE / reflectivity curves; f95 8.96 vs 8.42 Hz for the trapezoid)
  and the star temperature moves f90 / f99 by < 1 % and f95 by up to 40 % between
  5000 and 30 000 K only because r★ = 0.03 mas sits on the hump transition.
- **Sampling**: all stars with r★ ≥ 0.01 mas have f99 ≤ 34 Hz at 25 km/s, so 70 Hz is
  enough (`config.lightcurve_sampling`); searches binned to 5 ms (200 Hz) lose nothing.
  For noise numbers use the band 0.2–25 Hz (`analysis.bands.SIGNAL_BAND_HZ`).

---

## Searching an observation for occultations

`search_observation` runs the blind coincidence search on one recording:

```python
from kbo_occultation import load_observation_log
from kbo_occultation.config import BandpassConfig, GridConfig, NumericalConfig, StarConfig
from kbo_occultation.analysis.search import search_observation

log = load_observation_log()
obs = log.get("2025-12-16", "HD-17316")
dc = log.nights["2025-12-16"].dc_report

result = search_observation(
    obs.path,                               # .npz cache (a raw .bin also works)
    dc_reports={"m1": dc, "m2": dc, "lst1": dc},   # only for the DC-based corrections
    template_radii_m=(100.0, 300.0, 500.0),
    distance_au=43.0,
    star=StarConfig(temperature_K=30000, angular_radius_mas=0.03),
    band=BandpassConfig(300.0, 650.0, 25),
    grid=GridConfig(x_max_m=4000, n_x=800),
    numerics=NumericalConfig(n_r_grid=3000),
    correction="highpass",                  # raw | highpass | dc_detrend | dc_despike
    highpass_cutoff_hz=3.0,
    min_ntel=2,                             # telescopes required in coincidence
    n_background_shifts=100,                # time shifts for the accidental rate
)
for ev in result.coincidences:
    print(f"t={ev.time_s:.3f}s  n_tel={ev.n_tel}  combined SNR={ev.combined_snr:.1f}")
```

The pipeline stages are:
1. Load all three channels (variance flux proxy) at a fixed cadence.
2. Apply the noise correction per telescope.
3. Simulate a template dip for each KBO radius.
4. Slide each template over the light curve (matched filter) to get SNR(t).
5. Apply a χ² shape veto.
6. Keep events seen by at least `min_ntel` telescopes within the light-travel +
   cabling tolerance.
7. Estimate the accidental-coincidence rate from many time shifts.

On a benchmark star, where no KBOs are expected, any surviving coincidence is an
accidental, so these runs measure the pipeline's false-positive rate.
See [examples/real_data/search_observation.py](examples/real_data/search_observation.py).

---

### How the SNR is normalised

`sliding_matched_filter_snr(..., normalization=)` has three modes:

- `"empirical"` (default): divide the correlation by the rolling robust scatter of the
  filter output itself. On red noise the SNR then has σ ≈ 1 by construction; it is
  insensitive to how `compute_baseline_and_sigma` was run.
- `"local_block"`: the old white-noise normalisation (per-block MAD of the residual).
  It only sees noise above ~1/(block duration). On 1/f-like noise the SNR of pure
  noise has σ = 3.6 (exact median) to 5.3 (`decimate="auto"`, whose σ reads 30 % lower).
- `"global"`: one robust σ for the whole series.

`shape_veto_chi2(weighting="empirical")` weights the residuals with the same sigma, and
`calibrate_chi2_threshold` sets `max_chi2_reduced` from off-source lags. The injection
Monte Carlo normalises each short trial slice with the sigma measured once on the whole
uninjected series (measured on the slice it scatters ±10 % trial to trial).
Caveat: on a ~6 min run with noise coherent over > 10 s the empirical sigma is still
uncertain; SNR traces of real runs have σ ≈ 1.1–1.6 and MAD-σ ≈ 1.03
(`calibration.csv` of `kbo-night mimic`).

---

## How often does noise look like an occultation?

`kbo-night mimic 2025-12-16` (or `kbo mimic --night 2025-12-16`; none = all nights)
searches every logged star run with three corrections (raw, `dc_detrend`,
`dc_despike`), keeps all candidates above SNR 3 and writes `outputs/mimic/<nights>/`:

| File | Content |
|---|---|
| `candidates.csv` | per candidate: t, channel, SNR, χ²/dof, best R, depth, plus ±1 s features: DC slope and excursion flag of the target pixel, DC excursions of the background pixels, starguider deviation, wind (interpolated: those reports come every 10–16 s) |
| `events.csv` | coincidences (≥ 2 telescopes) of the real data, of the surrogates and the time-slide background |
| `far.csv`, `far.png` | false-alarm rate (events/hour) above each SNR: real vs Gaussian surrogates vs time slides, single-telescope and coincident |
| `tails.csv` | 3σ / 5σ / 7σ dip and bump counts, real vs surrogate vs white Gaussian, with the fraction overlapping a DC excursion |
| `flag_enrichment.csv` | fraction of real candidates carrying each auxiliary flag vs the fraction of random times carrying it |
| `calibration.csv` | σ of the SNR trace (should be ≈ 1) |

The surrogates (`analysis.surrogates`) are phase-randomised copies of each corrected
run: same spectrum, same gaps, no transients. A real rate above the surrogate rate
points to non-Gaussian mimics; a real coincident rate above the time-slide rate to
noise common to the telescopes.

First result, all three logged nights (good-time windows only: 0.73 telescope-hours,
0.21 h with ≥ 2 telescopes), events/hour at SNR ≥ 5, single telescope:

| correction | real | Gaussian surrogates |
|---|---|---|
| raw | 34 ± 7 | 15.6 ± 1.0 |
| dc_detrend | 26 ± 6 | 13.8 ± 1.0 |
| dc_despike | 23 ± 6 | 14.8 ± 1.0 |

The real rate is 1.6–2.2 × the Gaussian one. At SNR ≥ 5, 36 % of the raw candidates sit on a
DC excursion of the target pixel, against 1.4 % of random times (enrichment 27 ×);
background-pixel excursions and wind do not enrich them, starguider deviation does (3 ×).
`dc_despike` removes part of that excess but not all of it. Coincident events: none in
the real data (0.21 h), against 0.6–1.5 / h from time slides, so the coincidence
requirement works but the live time is too short to measure a coincident rate.
These numbers use 20 surrogates and 200 shifts per run and a 0.03 mas, 40 AU
template family; they will change as nights are added.

---

## Which noise removal costs the least detection power?

Figure of merit: efficiency ε(R, mag) at a fixed **real-data** false-alarm rate
([analysis/noise_removal.py](kbo_occultation/analysis/noise_removal.py),
[examples/studies/noise_removal_comparison.py](examples/studies/noise_removal_comparison.py),
tables in `outputs/noise_removal/`). For each method the SNR threshold is set so that the real
single-telescope FAR is 1 per telescope-hour, then star-only occultations (the sky background
and pedestal are not dimmed) are injected into the raw run, corrected, and searched with the
template as the method reshapes it. Methods: Butterworth cutoff scan (0.3–3 Hz), PSD-whitened
matched filter, DC detrend window (2, 5, 10 s), DC FFT swap (crossover ≤ 0.8 × DC Nyquist:
0.3 Hz MAGIC, 0.15 Hz LST-1; raised-cosine transition; tapered record ends), and the DC
despike used as a veto. Cross-telescope common-mode subtraction is deliberately not offered:
a real event is itself common-mode within ≤ 8 ms.

Efficiency averaged over the 3 logged nights (0.73 telescope-h, 20 injections per run, channel and radius):

| method | R = 50 m | 100 m | 300 m | 1000 m |
|---|---|---|---|---|
| whitened | 0.00 | 0.41 | 1.00 | 1.00 |
| butter 1 Hz | 0.00 | 0.27 | 1.00 | 1.00 |
| butter 3 Hz | 0.00 | 0.03 | 1.00 | 1.00 |
| dc_swap | 0.00 | 0.00 | 0.86 | 1.00 |
| dc_detrend 5 s | 0.00 | 0.00 | 0.71 | 1.00 |
| raw / dc_despike veto | 0.00 | 0.00 | 0.42 | 0.99 |
| dc_detrend 10 s | 0.00 | 0.00 | 0.30 | 1.00 |

- A 3 Hz high-pass is too aggressive: it keeps only 64 % of the white-noise SNR of a 100 m event
  (R = 1000 m: 54 %), against 89 % at 1 Hz (`check_highpass_signal_safety.py`, template embedded
  in a 240 s record). The pipeline default `config.HIGHPASS_CUTOFF_HZ` is still 3 Hz; this table
  favours ≈ 1–2 Hz or whitening.
- The old DC swap crossed over at 1 Hz, above the 0.5 Hz DC Nyquist; it is now 0.3 Hz.
- Caveat: the thresholds (10–24) are set by the largest few of ≈ 0.7 telescope-hours of real
  candidates, so they are noisy; add nights and re-run before changing a default.

---

## Injection / recovery and upper limits

Detection efficiency is measured by injecting simulated occultations into real
light curves and running the same search on them:

- `run_injection_monte_carlo(lc, radii_m, n_trials_per_radius, …, method_specs)`
  works on a single telescope. It compares noise corrections (`method_specs`:
  name → `LightCurve -> LightCurve`).
  [real_data/injection_montecarlo.py](examples/real_data/injection_montecarlo.py),
  [real_data/injection_snr_comparison.py](examples/real_data/injection_snr_comparison.py).
- `run_array_injection_monte_carlo(lcs, ORM_ARRAY, …)` works on the full array with
  coincidence. `efficiency_curve(trials, method)` turns its trials into
  efficiency vs radius.
- `SurveyExposure(star_name, live_time_s, shadow_velocity_mps, distance_au, efficiency)`
  per observed star, then `cumulative_density_upper_limit(exposures, radii_m)` gives
  the TAOS-style upper limit on the cumulative KBO surface density.
  [real_data/array_mc_benchmark.py](examples/real_data/array_mc_benchmark.py) runs
  the whole chain.

All real-data examples choose their input with `NIGHT, OBSERVATION = …` at the top
and `load_observation_log().get(NIGHT, OBSERVATION)`. They write to `outputs/`.

---

## Configuration reference

Configuration is code-based (dataclasses in [kbo_occultation/config.py](kbo_occultation/config.py)):

| Dataclass | Key fields | Purpose |
|-----------|------------|---------|
| `KBOConfig` | `radius_m`, `distance_au`, `impact_parameter_m=0.0` | The occulting object |
| `StarConfig` | `temperature_K`, `angular_radius_mas` | Background star (Planck spectrum + disk size) |
| `BandpassConfig` | `lam_min_nm`, `lam_max_nm`, `n_lambda` | Optical band + spectral sampling |
| `GridConfig` | `x_max_m`, `n_x` | Shadow-plane spatial grid |
| `NumericalConfig` | `n_int=40`, `n_r_grid=3000`, `n_star_side=32` | Numerical resolution |
| `TelescopeConfig` / `ArrayConfig` | positions (ENU), channel ↔ telescope ↔ DC key | Array geometry & coincidence |

`ORM_ARRAY` is the default array: MAGIC-1 (ch C, `m1`), MAGIC-2 (ch A, `m2`), LST-1
(ch B, `lst1`). The per-telescope physics constants used by `Instrument` (mirror
area, QE, reflectivity, PMT excess noise, …) are in
[MAGIC-LST.yml](kbo_occultation/data/iact_reference_values/MAGIC-LST.yml) and
`simulation/instruments.py`.

### Star photons and magnitude system

`Instrument.num_photons_from_mag` computes the star's photo-electrons as

    N = t · A · 10 · ∫ R(λ) Φ★(λ) dλ · 10^(−0.4 B)      (λ in nm, ×10 for nm→Å, A in cm²)
    Φ★ = Φ_Vega(λ) · 10^(0.4·0.03)

with R the full response (atmosphere × mirror × QE × funnel × camera × optical filter, if any) and Φ_Vega the
CALSPEC spectrum `alpha_lyr_stis_011.fits` in photons. Magnitudes are **Johnson B, Vega system** (Vega has B = 0.03),
whatever the optical filter; the no-filter and SII responses both lie inside Johnson B. The AB system (flat f_ν) is not
used: B_AB − B_Vega ≈ −0.15 mag and the reference shape differs (AB/Vega photons = 1.13 with no filter, 0.84 with SII).
Pass `temperature_K` (to `num_photons_from_mag`, `noise_budget`, ...) for a blackbody star normalised to the same
Bessell B photon count instead of Vega colours.

B = 0 star at airmass 1 (La Palma), effective Φ0 in ph cm⁻² s⁻¹ Å⁻¹:

| Telescope | Filter | photo-electrons / s | Φ0 |
|-----------|--------|--------------------:|---:|
| MAGIC | none | 6.83·10¹¹ | 1144 |
| MAGIC | `MAGIC_SII` | 1.16·10¹¹ | 1552 |
| LST | none | 1.73·10¹² | 1144 |
| LST | `MAGIC_SII` | 2.95·10¹¹ | 1552 |

Bundled data (`kbo_occultation/data/`): `optical_filters/` (MAGIC QE / SII / narrow,
Gaia, SDSS *ugriz*, Johnson `bessell_B`), `Vega_spectrum.txt` and `alpha_lyr_*.fits`
(`alpha_lyr_stis_011.fits` is the CALSPEC Vega spectrum used for the star photon counts), and
`iact_reference_values/` (atmospheric transmission, NSB spectra, telescope
reflectivity / QE curves).

---

## Other examples

[examples/studies/2025_12_16/](examples/studies/2025_12_16/) keeps the noise studies
of the first night (MAGIC-2 only) that led to the current corrections: real-noise
PSDs and noise vs bin size, high-pass filtering and its signal safety, baseline
fits, DC combination and DC-based noise removal. They are kept for reference; only
their imports are maintained.

---

## Running the tests

```bash
pytest            # runs tests/unit (see pyproject.toml)
```

---

## Known limitations & TODOs

See [PATCH_NOTES.md](PATCH_NOTES.md) for the history of fixes.

- **Telescope positions** in `ORM_ARRAY` are approximate placeholders; replace them
  with surveyed ENU coordinates when available (`TODO(user)` in `config.py`).
- **MAGIC optics** in `Instrument` reuse the LST QE / reflectivity curves scaled to
  the MAGIC mirror area, so MAGIC and LST differ only in collection area.
- **Per-frame noise floor**: `noise_budget` includes the variance-estimator floor
  √(2k/N) (k = 2, N = 2¹⁸ → 0.39 % per frame, MAGIC only, averaging down with the
  number of frames per bin) and an optional red-noise term `red_noise_rms` (default 0,
  not reduced by binning or by combining telescopes). `scint_coeff` defaults to 1.5.
- **Legacy functions** `compute_lightcurve_old` and `apply_stellar_disk` in
  `simulation/engine.py` are kept only for comparison and are unused/buggy. Use
  `compute_lightcurve` / `OccultationEngine`.

---

## Citation & license

If you use this software in a publication, please cite the authors
(E. do Souto Espiñeira, T. Hassan, V. Pascual) and the underlying occultation-survey
formalism referenced in `analysis/upper_limits.py` (Nihei et al. 2007; Zhang et al.
2008, 2013).

Released under the **MIT License**; see [LICENSE](LICENSE).
