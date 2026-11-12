# Minimum detectable KBO size: status and plan for steps 2 and 3

Goal: the smallest KBO that MAGIC-1 + MAGIC-2 + LST-1 can detect.
**Step 1 is done** (`kbo_occultation/analysis/min_size.py`, `tests/unit/test_min_size.py`).
Steps 2 and 3 are for a later session.

## Decisions (confirmed 2026-11-12)

- R_min = radius where the detection probability is **50 %**, with **≥ 2 of 3** telescopes in
  coincidence, at a coincidence FAR of **1 per hour** (0.1 and 0.01 per hour reported alongside).
- Main answer from the **measured** noise; the theory-only `noise_budget` is the "ideal system" reference.
- Targets: generic grid **B = 7 to 9**, `r★ = 0.03 mas · 10^(−0.2 (B − 7))`, Δ = 39 AU (a = 40 AU),
  opposition velocity, b = 0; **LST-1 infinity focus** scenario (star signal ×2); a loader for the
  user's **catalogue of suitable stars** (not the stars observed so far).

## What exists (inventory)

| Need | Where |
|---|---|
| Fresnel pattern, polychromatic + stellar disk | `simulation/physics.py`, `simulation/engine.py` (`OccultationEngine`, caches I(r) per R, Δ) |
| Instrument response, photons, white noise | `simulation/instruments.py` (`Instrument.response`, `noise_budget`) |
| Injection (star-only dilution) | `analysis/injection.py` (`inject_occultation(star_fraction=)`) |
| Blind search, empirical SNR normalisation, χ² veto | `analysis/matched_filter.py` |
| Coincidence, time-slide accidentals | `analysis/coincidence.py` |
| Real-data FAR, mimic candidates | `analysis/mimic.py` |
| Single-telescope ε(R) at fixed real FAR, per correction | `analysis/noise_removal.py` (`run_comparison`), results in `outputs/noise_removal/` |
| Array injection + coincidence (old; fixed threshold, dims the NSB too) | `analysis/injection_batch.py` (`run_array_injection_monte_carlo`) |
| Rates, observing time, density limits | `analysis/rates.py`, `analysis/upper_limits.py` |
| **Step 1 predictor** | `analysis/min_size.py` |

## Step 1: what was built

Frequency-domain prediction, no injections. Expected matched-filter SNR of the template in the
measured noise PSD (`spectra.matched_filter_sigma`) after a linear correction (`raw`, `whitened`,
`optimal`, `butter_<fc>`), in star units; per-run measured excess PSD (measured minus the model's
white level) added to the model white level of the target; thresholds from the real-data
candidate rates (`rate_curves`, `single_threshold`, `coincidence_threshold`, exponential tail
beyond the empirical range); Poisson-binomial ≥ k-of-n probability (`detection_probability`);
`r50`; `min_size_table`, `summarize`. Targets: `TargetStar`, `generic_targets`,
`load_target_catalogue`. Cache: `save_profiles` / `load_profiles`.

Typical use:

```python
from kbo_occultation.observations.log import load_observation_log
from kbo_occultation.analysis import min_size as ms

log = load_observation_log()
cfg = ms.MinSizeConfig()                                   # whitened, FAR 1/h, 2 of 3
profiles = ms.measure_noise_profiles(log, ["2026-08-19"], cfg)      # slow: loads the light curves
curves = ms.rate_curves("outputs/noise_removal/real_candidates.csv", log, ms.csv_method_name(cfg.correction))
table = ms.min_size_table(ms.generic_targets(), profiles, cfg, curves, focus_gains=(1.0, 2.0))
ms.summarize(table)
```

### First results (2026-08-19 profiles, whitened, FAR 1/h, not calibrated)

Median R50 [m] over the 4 runs of the night (range in the table), B = 7 / 8 / 9:

| configuration | B = 7 | B = 8 | B = 9 |
|---|---|---|---|
| MAGIC-2 | 99 | 96 | 92 |
| MAGIC-1 | 122 | 111 | 104 |
| LST-1 | 195 | 153 | 124 |
| 2 of 3 | 106 | 90 | 82 |
| 3 of 3 | 102 | 84 | 75 |
| theory-only, 2 of 3 | 42 | 43 | 52 |

Single-telescope thresholds at 1/h: A 9.3, B 11.2, C 8.3; 2-of-3 coincidence threshold 5.4 per
telescope (all extrapolated beyond the measured tail). In the measured case R50 falls with B because the excess-dominated noise does not grow with
magnitude while the star shrinks (sharper pattern, larger depth: to be confirmed); in the
theory case the photon loss wins and R50 rises.

### Known limitations found while verifying (what step 2 must calibrate)

1. **The prediction is pessimistic against injections.** Compared with
   `outputs/noise_removal/efficiency.csv` (same star, 0.03 mas, 40 AU, threshold from
   `thresholds.csv`): 50 m → 0 (agrees), ≥ 300 m → 1 (agrees except LST-1 at 8.6 mag, predicted 0,
   measured 1), but at 100 m the MAGIC channels are recovered 80-100 % where the prediction gives
   SNR 6-7 against a threshold of 10.9, i.e. the real SNR is ~1.5× higher (whitened); for `raw` the
   factor is 1.5-4. Likely cause: the real search normalises with a robust (MAD) scatter, which on
   intermittent red noise is below the rms that a Welch PSD gives; secondary: best-of-4 templates, the
   320-450 nm trapezoid band of the injection templates vs the instrument response used here, whitening measured at full cadence. `MinSizeConfig.snr_scale`
   multiplies every expected SNR for this; it is 1 until calibrated.
2. **LST focus has almost no effect in the measured case**, because the excess noise is held fixed in
   star units (`excess_scaling="star"`), so doubling the star light only lowers the white part.
   `excess_scaling="total"` (excess fixed in total-signal units) gives the optimistic bound. Which is
   right needs the HD-45563 focus data (or a focused run) to say how the excess scales.
3. The excess noise includes any NSB the model misses: the measured star fraction of LST-1 is
   0.27-0.68 against 0.65-0.9 in the model (HD-170263: 0.27 vs 0.66). It is carried to other
   magnitudes as a constant in star units.
4. Thresholds beyond the measured tail are extrapolated (0.7 telescope-hours); the coincidence
   threshold (5.4) is below the single-telescope threshold, so the 2-of-3 result is only as good as
   the exponential tail of three short candidate lists.
5. The detrend of the search (20 template durations) is not modelled; the DC-based corrections
   (`dc_swap`, `dc_detrend`) are not linear time-invariant and are not available.

## Step 2: empirical array efficiency, to validate and calibrate step 1

1. In `analysis/noise_removal.py`, add `injection_efficiency_array(...)`: inject the same star-only
   event (`star_fraction_for(obs, telescope)` per telescope) into channels A, B, C with the shadow
   offsets (`direction_enu` from `kinematics`; isotropic if unknown, as `injection_batch` does), correct
   each channel with the method (reshaped templates, as `reshaped_templates`), search each with the
   per-method threshold, then `coincidence.match_coincidences(min_ntel=2)` within
   `coincidence_tolerance_s`. Keep the affine segment speed-up.
2. Thresholds: per-telescope threshold from `min_size.coincidence_threshold` (common value), and check
   that the coincidence FAR of the real data with that threshold is the target, with
   `coincidence.accidental_coincidence_rate` (time slides) on the 3-channel nights.
3. Radii: about 12 values from 50 to 400 m; ≥ 50 trials per run, channel and radius; nights
   2026-08-19 and 2026-03-08. Regenerate `outputs/noise_removal/real_candidates.csv` with all nights
   first (the current one lacks 2026-03-08; `rate_curves` skips nights not in the CSV).
4. Compare measured ε(R) with `min_size_table`'s probability per run; fit `snr_scale` (and
   `snr_scatter`) per correction; decide `excess_scaling`; report the residual per run.
5. Do not change `injection_batch.py` (old weaknesses: dims the NSB, fixed threshold, local-block χ²).

## Step 3: study script, plots, rates

- `examples/studies/min_detectable_size.py`:
  - R50 vs B (7 to 9) for MAGIC alone, LST alone, 2 of 3, 3 of 3; measured vs theory; LST focus on/off
    (both `excess_scaling` values as a band); FAR 1, 0.1, 0.01 per hour;
  - `--catalogue file.csv --columns name=...,mag_B=...` for the user's catalogue (`load_target_catalogue`);
  - writes `outputs/min_size/` (table, summary, figures, the profiles `.npz`).
- Convert R50 and ε(R) into observing time with `rates.event_rate` / `upper_limits.SurveyExposure`,
  replacing the detection model of `tests/observing_time_vs_kbo_size.py` (shot-noise-only SNR at
  100 ms, a threshold on |I − 1|, three telescopes as √3).
- Extras: b-averaged efficiency inside H(r); `n_star_side` convergence (64 here; 128 was needed
  for converged f99, not for SNR); a template mismatch test (injected R ≠ template R, the search
  uses several radii); README section ("Minimum detectable size") and a notebook cell.
- If a catalogue star has a known spectral type, `TargetStar.temperature_K` sets the stellar colour
  in both the template and the photon rate (`model_budget` takes it from the star).
