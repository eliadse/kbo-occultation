# Patch notes

Files here mirror your repo's layout — drop them into `kbo_occultation/` and
`examples/` in place of the existing ones (or diff them against your
current versions first). Tested against a fresh clone of `main`.

## Bugs fixed

1. **Silent wrong results when reusing an `OccultationEngine` across impact
   parameters or star sizes** (`simulation.py`). `_get_fresnel_profiles`
   cached radial profiles by `(R_m, D_m, len(r_grid_m))`. Since
   `numerics.n_r_grid` is constant, two calls with *different*
   `impact_parameter_m` or star size (which change the needed `r_max`, not
   the grid length) hit the same cache key and got back a profile computed
   for the wrong radial range — no error, just a wrong light curve. Fixed
   by giving each engine a **fixed radial grid**, sized once at
   construction (auto-extended with a warning if a call needs more range),
   and caching by `(R_m, D_m)` alone, which is now safe.

2. **`NameError` in `OccultationEngine.compute()`'s point-source branch**
   (`simulation.py`): referenced an undefined `r_grid_m` instead of
   `r_obs`. Any call with `star.angular_radius_mas < 0.0001` crashed.

3. **`Instrument` class couldn't be instantiated** (`instruments.py`):
   used `ascii.read`, `interpolate.interp1d`, `u.deg`, and
   `reference_files_path` without importing `astropy.io.ascii`,
   `scipy.interpolate`, `astropy.units`, or defining that path. Added the
   imports and pointed `reference_files_path` at `PACKAGE_DATA`. **This
   only fixes the `NameError`/`AttributeError`** — the actual calibration
   files it reads (mirror reflectivity, QE curves, camera transmission,
   atmospheric transmission, NSB spectra — everything in
   `mirror_reflectivity_per_telescope`, `qes_per_telescope`,
   `camera_transmission_per_telescope`, plus the `atm_trans_*` and
   `Spectra_NSB_*` files) aren't in `kbo_occultation/data/` yet (only the
   two MAGIC filter files are). `Instrument()` will now fail with a clear
   `FileNotFoundError` naming the missing file instead of a confusing
   `AttributeError`, but it needs those files added before it actually
   works. I didn't fabricate calibration numbers — you'll know best which
   files are the right ones. Added `astropy` to `pyproject.toml`
   dependencies since it was used but never declared.

4. **Hard-coded absolute weight cutoff (`w < 1e-12`)** in
   `simulate_poly_point` and `OccultationEngine.compute`'s point-source
   branch: for narrow filters every individual wavelength weight can be
   well below `1e-12` even though it's a large *fraction* of the total,
   which would silently drop the whole bandpass. Replaced with a relative
   tolerance (`weight_rel_tol`, default `1e-6` of the max weight in the
   band) — this was already flagged as a to-do.

## New: `sweep.py`

`run_parameter_sweep(radii_m, distances_au, impact_parameters_m,
star_angular_radii_mas, star_temperature_K, bands, grid, numerics, ...)`
groups by bandpass and reuses one `OccultationEngine` per group across
every KBO radius/distance/impact-parameter/star-size combination, so the
Fresnel-profile cache (fixed by #1 above) actually pays off. Returns a
list of `SweepPoint` (params + depth/min_intensity + optionally the full
curve). `to_records()`/`to_dataframe()` give you a flat table for
plotting/filtering.

Only star *temperature* is a single fixed value per call — if you want to
sweep spectral type too, loop `run_parameter_sweep` over temperatures (it
already can't be cached across temperatures anyway, since that changes the
spectral weights).

## New: `detectability.py`

Noise-model-agnostic helpers that take a light curve and a per-sample
sigma (fractional flux) and return an SNR: `peak_snr` (deepest single
sample) and `matched_filter_snr` (whole-dip, better for wider/slower
events). Plus `spatial_to_time`/`resample_to_cadence` to convert the
spatial axis to a real sampling cadence via shadow velocity, and
`event_duration` / `is_detectable`. Deliberately doesn't compute sigma
itself — plug in the three-component noise model once it's ready, or a
flat placeholder for now (see `examples/parameter_sweep_example.py`).

## Not touched, but worth knowing about

- `simulation.py`'s `compute_lightcurve_old` and `apply_stellar_disk`
  (1D chord version) are still there, unused by the new code, with their
  own pre-existing issues (e.g. `compute_lightcurve_old` references
  `intensity_total` which is only defined inside an `if` with no `else`).
  Left alone since you may still want them for comparison, per your notes
  on the 1D-vs-2D convolution comparison.
- `photometry.py` has a `print("t0 = ", t_0)` typo (undefined `t_0`,
  should be `t0`) in `LightCurve.from_stat_binary`, and an
  `average_chunks(x, y, n_samples)` variant that's accidentally nested
  dead code inside `plot_lightcurves` (indentation puts it after a
  `return`, so it's unreachable — the module-level 2-arg
  `average_chunks(x, n)` right after it is the one that actually runs).
  Left both alone since they're outside what you asked about, but flagging
  in case they bite you later.
