# kbo-occultation

Simulation and detection of stellar occultation diffraction patterns by sub-km
Kuiper Belt Objects (KBOs), for fast photometry with the **MAGIC + LST-1** array at
the Roque de los Muchachos Observatory (ORM).

**Authors:** E. do Souto Espiñeira, T. Hassan, V. Pascual

---

## Overview

When a small KBO passes in front of a background star, it casts a Fresnel
diffraction shadow that sweeps across the Earth in milliseconds. This package lets
you:

- **Simulate** the diffraction light curve for a given KBO, star, and optical band —
  no telescope data required.
- **Search** your own fast-photometry recordings for these events with a blind
  matched-filter pipeline, and combine the three telescopes (MAGIC-1, MAGIC-2, LST-1)
  in coincidence to suppress false positives.
- **Characterise** sensitivity and set upper limits on the cumulative KBO surface
  density via signal injection / recovery Monte Carlo.

The physics core uses the Fresnel–Lommel (Hankel) formulation, polychromatic
integration over a Planck spectrum weighted by the filter response, and finite
stellar-disk convolution.

### Features

- Fresnel diffraction (Lommel–Hankel formulation)
- Polychromatic integration (Planck + filter response)
- Configurable astrophysical parameters (KBO size/distance, star temperature/size, band)
- Blind matched-filter search with a χ² shape veto
- Signal injection / recovery Monte Carlo (single telescope and full array)
- 3-telescope coincidence matching (MAGIC-1 / MAGIC-2 / LST-1)
- Upper limits on the cumulative KBO surface density

---

## Installation

Requires Python 3.9+ and the scientific-Python stack (numpy, scipy, astropy,
matplotlib, pandas — installed automatically).

```bash
git clone <repo-url>
cd kbo_occultation
pip install -e .
```

To also install the test dependency (pytest):

```bash
pip install -e ".[test]"
```

`kbo-occultation` is a **library** — there is no command-line tool. You drive it from
Python, either interactively or through the scripts in [examples/](examples/).

---

## Quickstart: simulate a light curve (no data needed)

Everything you need to produce a diffraction light curve lives in a handful of
dataclasses plus `compute_lightcurve`:

```python
import matplotlib.pyplot as plt
from kbo_occultation import (
    compute_lightcurve, KBOConfig, StarConfig,
    BandpassConfig, GridConfig, NumericalConfig,
)

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

See [examples/basic_simulation.py](examples/basic_simulation.py) for a monochromatic /
polychromatic / finite-star comparison, and
[examples/parameter_sweep_example.py](examples/parameter_sweep_example.py) to scan many
KBO radii, distances, impact parameters, and star sizes efficiently with
`run_parameter_sweep`.

---

## Working with your own data

The package processes fast-photometry recordings from the MAGIC + LST-1 array. The
three DAQ channels map to telescopes as **C = MAGIC-1, A = MAGIC-2, B = LST-1**. The
per-sample **flux proxy is the variance** (`std²`) of the digitised samples in each
record, not a mean count.

The real-data workflow has three steps:

1. **Place your raw recordings.** Drop the raw stat-binary files into
   `kbo_occultation/data/observations/`, named `Spectrum*.bin`. (This directory is
   git-ignored — observation data is never committed.) Optional slow "DC" reports go
   under `kbo_occultation/data/observations/DCs/<date>/dc_report.pkl`.

2. **Pre-process into compact caches.** Run
   [examples/preprocess_observations.py](examples/preprocess_observations.py). Each
   ~140 MB `.bin` is parsed once, the per-channel raw and outlier-cleaned variance plus
   the start time are extracted, and a small `.npz` cache is written beside the source.
   Re-run this after new data lands or after changing the outlier parameters.

3. **Run the blind coincidence search.** See
   [examples/search_observation_example.py](examples/search_observation_example.py),
   which calls `search_observation()`:

   ```python
   from kbo_occultation import PACKAGE_DATA
   from kbo_occultation.config import BandpassConfig, GridConfig, NumericalConfig, StarConfig
   from kbo_occultation.search import search_observation

   result = search_observation(
       f"{PACKAGE_DATA}/observations/<your_file>.npz",   # accepts .npz or raw .bin
       template_radii_m=(100.0, 300.0, 500.0),
       distance_au=43.0,
       star=StarConfig(temperature_K=30000, angular_radius_mas=0.03),
       band=BandpassConfig(300.0, 650.0, 25),
       grid=GridConfig(x_max_m=4000, n_x=800),
       numerics=NumericalConfig(n_r_grid=3000),
       correction="highpass",       # raw | highpass | dc_detrend | dc_despike
       highpass_cutoff_hz=3.0,
       min_ntel=2,                  # telescopes required in coincidence
       n_background_shifts=100,     # time-shifts for the accidental-rate estimate
   )

   for ch, res in result.per_telescope.items():
       print(res.telescope, ch, len(res.candidates), "candidates")
   for ev in result.coincidences:
       print(f"t={ev.time_s:.3f}s  n_tel={ev.n_tel}  combined SNR={ev.combined_snr:.1f}")
   ```

`search_observation()` accepts either a `.npz` cache or a raw `.bin` transparently.
Running it on a **benchmark star** (where no KBOs are expected) is how the pipeline's
false-positive rate is measured — any surviving coincidence there is an accidental.

---

## Pipeline stages

`search_observation()` (in [kbo_occultation/search.py](kbo_occultation/search.py))
performs, for one observation:

1. **Load** all three channels (variance flux proxy) at a fixed cadence.
2. **Correct** the noise per telescope (`raw`, `highpass`, `dc_detrend`, or `dc_despike`).
3. **Build templates** — simulate the expected diffraction dip for each candidate KBO radius.
4. **Matched filter** — slide each template over the light curve to get an SNR series.
5. **χ² shape veto** — reject candidates whose shape doesn't match the template.
6. **Coincidence** — keep events seen by at least `min_ntel` telescopes within the
   light-travel + cabling tolerance.
7. **Accidental background** — repeat the coincidence step over many time shifts to
   estimate the false-positive rate.

The same building blocks are exposed individually — `run_injection_monte_carlo` /
`run_array_injection_monte_carlo`, `efficiency_curve`, and
`cumulative_density_upper_limit` — to turn a detection efficiency into a KBO
surface-density upper limit (TAOS-style formalism).

---

## Configuration reference

Configuration is code-based (dataclasses in [kbo_occultation/config.py](kbo_occultation/config.py)),
not YAML files passed on a command line:

| Dataclass | Key fields | Purpose |
|-----------|------------|---------|
| `KBOConfig` | `radius_m`, `distance_au`, `impact_parameter_m=0.0` | The occulting object |
| `StarConfig` | `temperature_K`, `angular_radius_mas` | Background star (Planck spectrum + disk size) |
| `BandpassConfig` | `lam_min_nm`, `lam_max_nm`, `n_lambda` | Optical band + spectral sampling |
| `GridConfig` | `x_max_m`, `n_x` | Shadow-plane spatial grid |
| `NumericalConfig` | `n_int=40`, `n_r_grid=3000`, `n_star_side=32` | Numerical resolution |
| `TelescopeConfig` / `ArrayConfig` | positions (ENU), channel↔telescope map | Array geometry & coincidence |

`ORM_ARRAY` (in `config.py`) is the default array: MAGIC-1 (ch C), MAGIC-2 (ch A),
LST-1 (ch B) at the ORM. Import and pass a different `ArrayConfig` to override it. The
richer per-telescope physics constants (mirror area, QE, optical efficiency, PMT
excess noise, background) live in
[kbo_occultation/data/iact_reference_values/MAGIC-LST.yml](kbo_occultation/data/iact_reference_values/MAGIC-LST.yml).

---

## Data & required files

**Bundled with the package** (`kbo_occultation/data/`):

- `optical_filters/` — MAGIC QE / SII / narrow-filter curves, plus Gaia and SDSS *ugriz* responses.
- `Vega_spectrum.txt`, `alpha_lyr_*.fits` — reference stellar spectra.
- `iact_reference_values/` — atmospheric transmission, NSB spectra, telescope
  reflectivity/QE curves, and `MAGIC-LST.yml` (array geometry + physics constants).

**Supplied by you** (not shipped, git-ignored):

- `data/observations/Spectrum*.bin` — raw stat-binary recordings. The format is defined
  in [kbo_occultation/io.py](kbo_occultation/io.py) `read_stat_binary_file()`: a record
  array of `mean/std/min/max` for channels A–D plus a `time_stamp`. Only `std` (→
  variance, the flux proxy) for channels A/B/C and the start time are used; timestamps
  are rebuilt at a fixed cadence.
- `data/observations/.../*.npz` — compact caches written by the pre-processing step
  (raw + cleaned per-channel variance, `t0`, `dt`, `n`, JSON metadata).
- Optional `data/observations/DCs/<date>/dc_report.pkl` — slow "DC" reports for the
  `dc_detrend` / `dc_despike` corrections.

---

## Examples index

All runnable scripts live in [examples/](examples/):

**Simulation**
- [basic_simulation.py](examples/basic_simulation.py) — monochromatic vs. polychromatic vs. finite-star light curves.
- [plot_example.py](examples/plot_example.py) — minimal single-curve plot.
- [parameter_sweep_example.py](examples/parameter_sweep_example.py) — scan a grid of parameters with `run_parameter_sweep`.
- [plot_diffraction_psd.py](examples/plot_diffraction_psd.py) / [plot_lightcurve_comparison.py](examples/plot_lightcurve_comparison.py) — diffraction power spectrum & comparisons.
- [magic_lst_filter_comparison.py](examples/magic_lst_filter_comparison.py) — compare the MAGIC and LST optical bands.

**Injection / recovery Monte Carlo**
- [injection_montecarlo_2025_12_16.py](examples/injection_montecarlo_2025_12_16.py) — single-telescope recovery efficiency and false-alarm study.
- [injection_snr_comparison_2025_12_16.py](examples/injection_snr_comparison_2025_12_16.py) — SNR comparison across methods.
- [array_mc_benchmark.py](examples/array_mc_benchmark.py) — full-array (coincidence) Monte Carlo.

**Noise studies**
- [inspect_real_noise_2025_12_16.py](examples/inspect_real_noise_2025_12_16.py), [highpass_filter_2025_12_16.py](examples/highpass_filter_2025_12_16.py), [check_highpass_signal_safety.py](examples/check_highpass_signal_safety.py), [baseline_fit_test_2025_12_16.py](examples/baseline_fit_test_2025_12_16.py), [dc_combine_2025_12_16.py](examples/dc_combine_2025_12_16.py), [dc_noise_removal_comparison_2025_12_16.py](examples/dc_noise_removal_comparison_2025_12_16.py).

**Real-data search**
- [preprocess_observations.py](examples/preprocess_observations.py) — build `.npz` caches from raw `.bin`.
- [search_observation_example.py](examples/search_observation_example.py) — end-to-end blind coincidence search.

---

## Running the tests

```bash
pip install -e ".[test]"
pytest
```

The default suite lives in [tests/unit/](tests/unit/) (configured via `testpaths` in
[pyproject.toml](pyproject.toml)). The other scripts under [tests/](tests/) are
longer-running physics / injection-recovery checks that produce diagnostic plots.

---

## Known limitations & TODOs

See [PATCH_NOTES.md](PATCH_NOTES.md) and [TODO.txt](TODO.txt) for details. In short:

- **Telescope positions** in `ORM_ARRAY` are approximate placeholders; replace with
  surveyed ENU coordinates when available (`TODO(user)` in `config.py`).
- **`instruments.Instrument`** needs several calibration files (reflectivity, QE, camera
  transmission, atmospheric transmission, NSB spectra) that are not all bundled yet; it
  raises a clear `FileNotFoundError` naming the missing file until they are added.
- **Legacy functions** `compute_lightcurve_old` and `apply_stellar_disk` in
  `simulation.py` are kept only for comparison and are unused/buggy — prefer
  `compute_lightcurve` / `OccultationEngine`.

---

## Citation & license

If you use this software in a publication, please cite the authors
(E. do Souto Espiñeira, T. Hassan, V. Pascual) and the underlying occultation-survey
formalism referenced in `upper_limits.py` (Nihei et al. 2007; Zhang et al. 2008, 2013).

Released under the **MIT License** — see [LICENSE](LICENSE).
