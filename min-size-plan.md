# Minimum detectable KBO size: status and plan for step 3

Goal: the smallest KBO that MAGIC-1 + MAGIC-2 + LST-1 can detect.
**Steps 1 and 2 are done** (`kbo_occultation/analysis/min_size.py`,
`examples/studies/min_size_calibration.py`, `tests/unit/test_min_size.py`, `tests/unit/test_noise_removal.py`).
Step 3 is for a later session.

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
cfg = ms.MinSizeConfig().calibrated()                      # whitened, FAR 1/h, 2 of 3, SNR_CALIBRATION (step 2)
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

### Known limitations found while verifying (step 2 outcomes in italics)

1. **The prediction is pessimistic against injections.** Compared with
   `outputs/noise_removal/efficiency.csv` (same star, 0.03 mas, 40 AU, threshold from
   `thresholds.csv`): 50 m → 0 (agrees), ≥ 300 m → 1 (agrees except LST-1 at 8.6 mag, predicted 0,
   measured 1), but at 100 m the MAGIC channels are recovered 80-100 % where the prediction gives
   SNR 6-7 against a threshold of 10.9, i.e. the real SNR is ~1.5× higher (whitened); for `raw` the
   factor is 1.5-4. Likely cause: the real search normalises with a robust (MAD) scatter, which on
   intermittent red noise is below the rms that a Welch PSD gives; secondary: best-of-4 templates, the
   320-450 nm trapezoid band of the injection templates vs the instrument response used here, whitening measured at full cadence. `MinSizeConfig.snr_scale`
   multiplies every expected SNR for this; it is 1 until calibrated.
   *Step 2: for `whitened` the gap was mostly an artefact of comparing with the single-telescope threshold (10.9). At the
   coincidence thresholds (about 5) the whitened prediction is already unbiased (scale 0.91). The pessimism is real for
   `raw` (1.54, LST-1 2.4) and on the faintest star (HD170263).*
2. **LST focus has almost no effect in the measured case**, because the excess noise is held fixed in
   star units (`excess_scaling="star"`), so doubling the star light only lowers the white part.
   `excess_scaling="total"` (excess fixed in total-signal units) gives the optimistic bound. Which is
   right needs the HD-45563 focus data (or a focused run) to say how the excess scales.
   *Step 2: the excess-vs-star-fraction regression gives slope -2.3 +- 1.0 (LST-1 -3.3 +- 0.8, MAGIC -0.8 +- 1.7; 0 =
   "star", -2 = "total"). It leans to "total" but does not decide; injections cannot decide it.*
3. The excess noise includes any NSB the model misses: the measured star fraction of LST-1 is
   0.27-0.68 against 0.65-0.9 in the model (HD-170263: 0.27 vs 0.66). It is carried to other
   magnitudes as a constant in star units.
4. Thresholds beyond the measured tail are extrapolated (0.7 telescope-hours); the coincidence
   threshold (5.4) is below the single-telescope threshold, so the 2-of-3 result is only as good as
   the exponential tail of three short candidate lists.
5. The detrend of the search (20 template durations) is not modelled; the DC-based corrections
   (`dc_swap`, `dc_detrend`) are not linear time-invariant and are not available.

## Step 2: array injections to validate and calibrate step 1 (done 2026-11-12)

Built:

- `noise_removal.py`: `shadow_direction_enu`, `shadow_offsets_s`, `injection_efficiency_array` (star-only injection into
  all channels with the shadow travel-time offsets, the reshaped bank search of the real data, any threshold and any k-of-n
  read off afterwards), `ComparisonConfig.edge_margin_s` / `peak_window_s`. `mimic.search_lightcurve(peak_near=)` returns the
  largest SNR near the injection, above the floor or not. `injection_efficiency`, `run_comparison` and `injection_batch.py`
  are untouched.
- `min_size.py`, section "calibration against injections": `calibrate_run`, `timeslide_far`, `thresholds_table`,
  `measured_efficiency`, `predicted_snr`, `predicted_efficiency`, `fit_snr_calibration`, `scale_by_deviance`,
  `excess_scaling_test`, `efficiency_r50`, `observed_target`, `calibration_runs`, `injection_patterns` (and `_engine`, split
  out of `template_family`), the fitted `SNR_CALIBRATION`, `MinSizeConfig.calibrated()`.
- `examples/studies/min_size_calibration.py` (resumable; writes `outputs/min_size/calibration/`):
  `python examples/studies/min_size_calibration.py 2026-03-08 2026-08-19 --trials 50 --corrections whitened raw butter_1`
  (about 30 min for the injections). Smoke test: one night, `--trials 5 --radii 50 100 300 --model-sf-trials 0`.

Differences from the plan as first written:

- `real_candidates.csv` did **not** need regenerating: it already has 2026-03-08 and its live time matches the log. The
  script checks that each run and method is present.
- `injection_efficiency` needs more than 120 s of good time (it re-corrects +-60 s around the event), so it gave NaN for
  5 of the 7 three-channel runs. The array version draws the event inside the good-time window and crops the search to it.
- The "model" star-fraction check uses `profile.star_fraction_model` (airmass, throughput and PSF of the run), not
  `star_fraction_for`.
- The scale is applied to measured noise only; the theory-only ideal system is not scaled.
- 2025-12-16 has one channel, so it is not in the array injections.

Setup of the results below: 2026-03-08 and 2026-08-19, 7 runs x 3 channels, 12 radii (50-400 m), 50 trials, measured star
fraction, the same thresholds for prediction and injections (FAR 1/h), whitened / raw / butter_1.

| correction | scale (measured / predicted SNR) | scatter | likelihood-fit scale | median R50 error, step 1 -> calibrated |
|---|---|---|---|---|
| whitened | 0.91 (per run 0.67-1.65) | 1.23 | 0.96 | +2 % -> +7 % (typical abs 9 % -> 7 %) |
| butter_1 | 1.08 (0.89-1.58) | 1.46 | 1.15 | +8 % -> +5 % |
| raw | 1.54 (0.87-6.5; LST-1 2.4) | 2.25 | 1.37 | +29 % -> +3 % (typical abs 8 %) |

Findings:

1. The step-1 pessimism is not general. Whitened is already right (scale 0.91, within the run-to-run spread), `butter_1` is
   close, `raw` is 1.5x low (LST-1 2.4x). The scale follows the full-run SNR-trace std (about 1.1) for whitened and
   `butter_1`, but not for `raw` on LST-1, so the robust-scatter explanation covers only part of it.
2. The faintest star (HD170263, B = 8.62) is the worst case: the whitened R50 prediction is 38-42 % too pessimistic
   (162 m predicted, 117 m measured), and LST-1 reaches 64 % at 400 m against 34 % (step 1) and 14 % (calibrated). Not understood yet; relevant
   because step 3 covers B = 7-9.
3. Thresholds: the time-slide threshold (2-of-3, 2000 shifts, 3 channel runs) agrees with the independent-stream one within
   about 0.1-0.6 SNR at 1/h and 0.1/h; at 0.01/h the analytic one is 0.5-1.3 higher (conservative).
4. Star fraction: the model fraction gives a whitened scale of 1.02 (LST-1 1.14) against 0.91 (0.93) for the measured one, so
   the choice moves LST-1 by about 20 %.
5. The peak SNR over the search window and the 4 templates is biased high at low SNR (at 50 m the median is 2.9 against 1.2
   predicted), so the scale fit uses cells with a median peak SNR of 5-50.

Calibrated step-1 table (whitened, FAR 1/h, 2 of 3, generic stars, median R50 [m] over the runs, B = 7 / 8 / 9):

| configuration | 2026-08-19 (4 runs) uncalibrated | calibrated | both nights (7 runs) calibrated |
|---|---|---|---|
| MAGIC-2 | 86 / 79 / 77 | 90 / 83 / 80 | 76 / 69 / 70 |
| MAGIC-1 | 122 / 111 / 103 | 128 / 116 / 109 | 105 / 94 / 91 |
| LST-1 | 189 / 146 / 127 | 202 / 157 / 136 | 170 / 128 / 115 |
| 2 of 3 | 103 / 89 / 80 | 107 / 92 / 84 | 81 / 73 / 70 |
| 3 of 3 | 104 / 84 / 75 | 109 / 89 / 80 | 93 / 72 / 65 |
| theory-only, 2 of 3 | 41 / 43 / 49 | (not scaled) | |

(`outputs/min_size/step1_table_calibrated.csv`; radii grid 40-500 m, 20 values.) LST focus x2 leaves these almost unchanged
in the measured case (`excess_scaling="star"`). 2026-03-08 is much quieter than 2026-08-19, which is why the two sets differ.

Open after step 2:

- `excess_scaling` is undecided (the "total" bound is the optimistic LST-focus case).
- The faint-star pessimism (finding 2), the coincidence thresholds beyond the measured tail, and the run-to-run spread of
  the scale (about +-30 %) are the systematic uncertainty on R50; a typical single-run error after calibration is 7-10 %,
  up to 40 %.
- Only two nights with 3 channels and short good-time windows (47-350 s).

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
