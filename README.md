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
| Add a new observation night | [`kbo-night scaffold` / `prepare`](#adding-a-new-night) |
| Measure the SNR of real observations and compare it with theory | [`kbo-night snr`](#measuring-the-snr-of-real-observations) |
| Simulate an occultation light curve or a theoretical SNR | [`compute_lightcurve`, `Instrument`](#simulating-light-curves-and-theoretical-snr) |
| Search a recording for occultations | [`search_observation`](#searching-an-observation-for-occultations) |
| Measure detection efficiency and set upper limits | [injection Monte Carlo](#injection--recovery-and-upper-limits) |

---

## Installation

Requires Python 3.9+ and the scientific-Python stack (numpy, scipy, astropy,
matplotlib, pandas, pyyaml — installed automatically).

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
  cli.py              the kbo-night command
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
| `observations.prepare` | `scaffold_night`, `prepare_night`: from raw `.bin` files to analysis-ready caches |
| `observations.photometry` | `LightCurve`, `.npz` pre-processing |
| `observations.io` | stat-binary file reader |
| `observations.filtering` | high-pass filters, outlier removal |
| `observations.dc_combine` | slow DC reports: loading, combining with / detrending / despiking the fast signal |
| `observations.cc_reports` | MAGIC CC `.rep` parsing, pointing / moon geometry, `conditions_table` |
| `analysis.snr` | `snr_report`: measured vs theoretical SNR per observation (what `kbo-night snr` runs) |
| `analysis.snr_measurement` | the individual SNR, NSB and noise estimators |
| `analysis.matched_filter`, `analysis.search`, `analysis.coincidence` | blind matched-filter search with χ² veto and multi-telescope coincidence |
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
   and the DC report in `DCs/<YYYY_MM_DD>/dc_report.pkl`. Put the CC `.rep` files anywhere
   under `cc_reports_root`.

2. **Generate the log entry.**

   ```bash
   kbo-night scaffold 2026-10-20
   ```

   This prints a ready-to-paste YAML block for the night. It covers every
   `Spectrum*.bin` whose filename timestamp falls between 12:00 UTC on that date and
   12:00 UTC the next day. Observation names, files, DAQ settings, dark runs (names
   containing "dark") and the DC report path are filled in from the filenames. Every
   other field is `null  # TODO(user)`. Paste the block under `nights:` in
   `observation_log.yaml`.

3. **Fill in what only you know.** Get coordinates, B/V magnitudes and spectral types
   from SIMBAD. Add the CC report paths, and remove channels that recorded no star
   light. The fields:

   | Field | Meaning |
   |---|---|
   | `channels` | DAQ channels that recorded star light, e.g. `[A]` if only MAGIC-2 was connected |
   | `dc_report`, `cc_reports` | DC pickle (relative to the data dir); CC `.rep` per telescope (`m1`, `m2`; relative to `cc_reports_root`) |
   | `dark_run` | name of the night's dark run / dark patch (NSB reference); `null` if none |
   | `faulty_mirrors` | `{m1: n, m2: n, lst1: n}` panels not correctly focused (AMC undefined + error); an observation's own `faulty_mirrors` overrides the night's |
   | `type` | `star` or `dark` |
   | `ra_deg`, `dec_deg`, `mag_B`, `mag_V`, `spectral_type` | from SIMBAD (ICRS degrees) |
   | `good_time` | `[start, stop]` seconds from the start of the run to use for the SNR; `null` = whole run |
   | `elog_notes`, `notes`, `manual`, `manual_conditions` | free text, and manual overrides of measured conditions |

   Zd, airmass and moon are **not** entered by hand. They are computed from the
   coordinates and the run times.

4. **Prepare the night.**

   ```bash
   kbo-night prepare 2026-10-20
   ```

   This builds the missing `.npz` caches, condenses the CC reports into
   `conditions/2026-10-20/`, prints the AMC unfocused-panel count of each run (copy
   these into `faulty_mirrors`), and lists everything still missing from the log.
   Existing caches are kept; `--force` rebuilds them, e.g. after changing the
   outlier-removal parameters. Re-run it until the missing list only contains things
   you know to be unavailable.

5. **Check and measure.** `kbo-night list` shows every night's status, and
   `kbo-night snr 2026-10-20` gives the SNR table (next section).

---

## Measuring the SNR of real observations

```bash
kbo-night snr 2026-10-20 --plot          # one night
kbo-night snr --plot                     # every night in the log
kbo-night snr 2026-10-20 --bins 10 50 200 --detrend-window 10
```

These commands write `outputs/snr/snr_<nights>.csv` with all columns, print a summary,
and with `--plot` save SNR-vs-magnitude figures per bin size. The same from Python:

```python
from kbo_occultation import load_observation_log, snr_report

log = load_observation_log()
results = snr_report(log, ["2026-10-20"], n_bins=(10, 200))
results.loc[("2026-10-20", "HD-12345", "MAGIC-1", 10), ["detr", "theo_obs_B", "ratio_detr_B"]]
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
| `*_ext` | × extinction factor 10^(0.4 k_B (X − 1)) | |

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
`kbo-night conditions NIGHT [--csv FILE]` prints the per-observation conditions on
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
magic.noise_budget(7.0, 5 * u.ms, airmass=1.2) # dict: sigma_shot, sigma_elec, sigma_scint, sigma_total, ...
band = BandpassConfig(*magic.response_band())  # band covering the instrument response, for compute_lightcurve
```

Examples in [examples/simulation/](examples/simulation/):
- [basic_simulation.py](examples/simulation/basic_simulation.py): monochromatic vs polychromatic vs finite star
- [plot_example.py](examples/simulation/plot_example.py): minimal single curve
- [parameter_sweep_example.py](examples/simulation/parameter_sweep_example.py): `run_parameter_sweep` and a flat-noise detectability cut
- [plot_lightcurve_comparison.py](examples/simulation/plot_lightcurve_comparison.py): light curves at the real cadence with instrument noise
- [plot_diffraction_psd.py](examples/simulation/plot_diffraction_psd.py): power spectrum of the diffraction pattern

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

Bundled data (`kbo_occultation/data/`): `optical_filters/` (MAGIC QE / SII / narrow,
Gaia, SDSS *ugriz*), `Vega_spectrum.txt` and `alpha_lyr_*.fits`, and
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
- **Per-frame noise floor**: the measured MAGIC per-frame scatter has a ≈0.4 %
  level-independent floor that is consistent with the variance-estimator error.
  `noise_budget` does not model it.
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
