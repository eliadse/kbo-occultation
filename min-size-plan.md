# Minimum detectable KBO size: status and plan for step 3

Goal: the smallest KBO that MAGIC-1 + MAGIC-2 + LST-1 can detect.
**Steps 1, 2, the follow-up "Phase A" and step 3 are done** (`kbo_occultation/analysis/min_size.py`,
`examples/studies/min_size_calibration.py`, `examples/studies/min_detectable_size.py`, `tests/unit/test_min_size.py`,
`tests/unit/test_rates.py`, `tests/unit/test_noise_removal.py`). What is left is under "Open".

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

## Phase A: open items of step 2 (done 2026-11-16)

Re-checking the step-2 outputs showed that HD170263 (the only long run, 347 s) has an SNR trace with std 1.12 but a
MAD-sigma of 0.89, i.e. intermittent noise: the search divides by a local robust scatter (20 template durations,
about 16 s) whereas the prediction used one Welch PSD for the whole run, so bursts raised the predicted noise everywhere.

- **Time-resolved noise (adopted, `MinSizeConfig.time_resolved = True`).** `measure_noise_profiles(blocks=True)` also
  gives every ~16 s block of a run its own median-Welch excess PSD (`NoiseProfile.S_excess_blocks`, `spectra.blockwise_welch`);
  `per_telescope_snr` returns one SNR row per block (the correction G stays that of the whole run) and `detection_probability`
  averages the probability over the blocks (an event falls at a uniform time). Checked against the cached step-2 injections,
  no new ones needed:

  | | stationary | time resolved |
  |---|---|---|
  | HD170263 2-of-3 whitened R50 error (calibrated) | +42 % | +14 % |
  | per-run scale spread, whitened / butter_1 / raw (std) | 0.22 / 0.20 / 1.43 | 0.14 / 0.11 / 0.23 |
  | whitened R50 error, mean abs (calibrated) | 9.0 % | 5.5 % |
  | whitened scale | 0.91 | 1.07 |

  The whitened scale is now about 1, i.e. the physical prediction is right without a fudge factor. Stationary results stay
  available (`--stationary`, `SNR_CALIBRATION_STATIONARY`; old outputs in `outputs/min_size/calibration_stationary_2026-11-12/`).
- **Calibration with 2025-12-16** (5 single-channel MAGIC-2 runs, `--min-channels 1`): 26 run-channels,
  `SNR_CALIBRATION = {whitened: (1.072, 1.174), raw: (0.655, 1.774), butter_1: (0.971, 1.221)}`. The whitened scale does
  not depend on magnitude (Spearman -0.00, p = 0.99) or on the measured star fraction (+0.17, p = 0.40): the faint-star
  pessimism of step 2 was HD170263-specific (burstiness), not a trend in B. Calibrated whitened 2-of-3 R50 error: median
  +1.7 %, mean abs 5.5 %, worst run +14 %.
- **Excess noise vs star fraction**, 26 profiles: slope -2.43 +- 0.89 (MAGIC -1.22 +- 1.47, LST-1 -3.33 +- 0.77), still
  leaning to "total" but not decisive; no focused LST-1 run exists in the log, so both stay as a band.
- Adding 2025-12-16 to the rate curves raises the MAGIC-2 threshold at 1/h from 6.6 to 9.3 and the 2-of-3 one from 5.17 to 5.39.
- `thr_at_floor` column (threshold equal to the 3.0 candidate floor, i.e. the FAR is only an upper bound): no row is flagged
  once 2025-12-16 is included (3-of-3 at 1/h: 3.30).
- Radius grid for step 3: 10-1000 m, 41 values (the theory-only R50 is no longer below the grid).

## Step 3: study script, plots, rates (done 2026-11-16)

Built: `examples/studies/min_detectable_size.py` (see its docstring; smoke test in it, full run ~1 h:
`python examples/studies/min_detectable_size.py --check-numerics`, writes `outputs/min_size/`: `r50_table.csv`,
`r50_summary.csv`, `thresholds.csv`, `efficiency.csv`, `detection_width.csv`, `detection_width_summary.csv`,
`observing_time.csv`, `observing_time_shot_noise_model.csv`, `numerics.csv`, `catalogue_r50.csv` with `--catalogue`, and
`r50_vs_mag.png`, `r50_lst_focus.png`, `efficiency_curves.png`, `observing_time.png`).

Library: `min_size.min_size_table` (`excess_scalings=`, `scale_band=`, `return_efficiency=`; `scale_quantiles`),
`expected_snr(signal=)` / `bank_snr` / `bank_residuals` (mismatched events and the real 4-template bank),
`snr_evaluator`, `array_threshold`, `parse_column_map`, `load_or_measure_profiles`, and the impact-parameter section
(`shadow_width_m`, `impact_nodes`, `width_signals`, `detection_width`); `rates.observing_time_per_size`,
`rates.detected_event_rate`; `upper_limits.EfficiencyCurve(clip=False)`; `tests/observing_time_vs_kbo_size.py --width-csv`
overlays the measured-noise widths on the old shot-noise plot.

### Results (whitened, calibrated, time-resolved; 7 three-channel runs; FAR 1/h; generic stars, 39 AU, b = 0)

R50 [m], median (range over the 7 runs) [calibration band, 16-84 % of the per-run scale]:

| configuration | B = 7 | B = 8 | B = 9 |
|---|---|---|---|
| MAGIC-2 | 93 (87-130) [87-96] | 80 (76-116) [75-83] | 78 (76-107) [74-81] |
| MAGIC-1 | 102 (90-171) [95-105] | 89 (80-154) [84-93] | 85 (78-140) [80-88] |
| LST-1 | 166 (108-365) [155-171] | 128 (88-287) [120-132] | 113 (80-244) [106-118] |
| **2 of 3** | **83 (72-135) [78-86]** | **72 (61-122) [68-75]** | **67 (58-111) [63-70]** |
| 3 of 3 | 91 (65-183) [86-94] | 74 (55-147) [70-77] | 64 (53-128) [60-66] |
| theory-only, 2 of 3 | 42 | 43 | 50 |

- FAR (2 of 3): 0.1/h gives 92 / 80 / 75 m, 0.01/h gives 101 / 87 / 83 m (thresholds per telescope 5.4 / 6.7 / 8.0, all extrapolated
  beyond the measured tail).
- R50 falls with B in the measured case because the excess-dominated noise does not grow with magnitude while the star shrinks
  (sharper pattern, deeper dip); in the theory-only case the photon loss wins and R50 rises (42 to 50 m).
- Excess scaling: with "total" (excess fixed in total-signal units) 2 of 3 gives 83 / 79 / 84 m, so at B = 9 the unknown scaling is
  worth +17 m (+25 %), more than the calibration band (+-5 %). It is the biggest systematic at faint magnitudes.
- LST-1 at infinity focus (star light x2) changes almost nothing: 2 of 3 unchanged with "star", 83 / 78 / 80 m with "total"; LST-1 alone
  166 / 127 / 111 ("star") or 162 / 131 / 124 m ("total"). The LST-1 contributes little to the 2-of-3 R50, which is set by the MAGICs.
- Numerics (`numerics.csv`): `n_r_grid` 4000 to 8000 changes the expected SNR by < 0.01 %, `n_star_side` 64 to 128 by <= 0.18 %.
- **The search bank is matched to a 0.03 mas star**, the prediction uses the target's own template. The expected SNR found by the
  4-template bank, against the perfect template, is 0.99 / 0.88 / 0.76 / 0.66 / 0.60 at B = 7 / 7.5 / 8 / 8.5 / 9 (the star is 0.012 mas at
  B = 9; the loss does not depend on R). Cutting the SNR by that factor (`outputs/min_size/bank_effect.csv`, 2 of 3, FAR 1/h) leaves R50 at
  84 / 83 / 83 / 84 / 87 m instead of 83 / 77 / 72 / 68 / 67 m: with the current bank the fall of R50 with B disappears.

Detection width W(R) = 2 integral P(R, b) db (2 of 3, measured, calibrated): at the b = 0 R50 the width is only ~10 % of H(R) (B = 7, R = 79 m:
P(b = 0) = 0.33, W = 0.5 km against H = 5.6 km); W / H reaches 0.5 at R about 160 m (B = 7), and exceeds 1 for B = 8-9 at large R (the bright
fringes outside the geometric shadow are detected; P at the outermost impact parameter is < 2e-11, so the integral is not cut short).
So R50 for a central chord flatters the cross-section.

Observing time per detected event, whole size distribution (Nihei q_s = 4, break 80 km, one star, 39 AU): 2 of 3 measured 273 / 183 / 146 h at
B = 7 / 8 / 9; theory-only 36 / 39 / 55 h; MAGIC-2 391 / 278 / 247 h; MAGIC-1 519 / 359 / 299 h; LST-1 2400 / 1352 / 787 h. For one size the minimum
for 2 of 3 is 572 / 387 / 335 h at D = 200 / 159 / 159 m (B = 7 / 8 / 9; theory-only about 70-115 h at D = 100-125 m). The old shot-noise model (3 MAGIC, 100 ms, threshold |I - 1| > 1/SNR) says a few hours
for D = 30-60 m: its widths (10-20 km) are larger than H(R) because the threshold is tiny, so it is very optimistic.

## Synthetic runs (2026-11-16)

`examples/studies/synthetic_min_size.py` (`analysis/synthetic.py`) runs the injection / blind search / threshold / efficiency / R50 chain on
generated Gaussian noise (white = noise model, coloured = + the measured excess of HD199102, 2026-08-19), stationary, independent channels,
stars of 0.03 mas · 10^(-0.2 (B - 7)), whitened. Medium run (`outputs/min_size/synthetic_quick`: B = 7 / 8 / 9, one 400 s run each, 10 trials x 8 radii
of 20-300 m, thresholds from the uninjected synthetic search itself, 0.33 h per channel and noise model, 2-of-3 thresholds extrapolated;
made with the **fixed 0.03 mas bank** of the older search, i.e. today's `--fixed-bank`):

- The search normalisation is right on this noise: full-run SNR-trace std 0.98-1.01 in every channel.
- Measured / predicted peak SNR against the **perfectly matched** template: 0.84 (white), 0.67-0.79 (coloured); against the **fixed 0.03 mas search bank**
  (`predicted_snr(bank=bank_residuals(cfg))`): 0.99-1.05 (white) and 0.99-1.00 (coloured), run-channel range 0.95-1.09. The likelihood fit on the 2-of-3
  detections gives 1.15 / 1.08 (white / coloured; few trials). So for stationary Gaussian noise the predictor, with the bank's template mismatch, is
  right to a few percent; the 0.91 of `SNR_CALIBRATION_STATIONARY` (0.03 mas stars, bank matched) and the 1.07 of the time-resolved one come from the real
  noise (intermittent, non-Gaussian), not from the search.
- 2-of-3 R50 [m] at B = 7 / 8 / 9, FAR 1/h: coloured measured 81 / 83 / 84 (bank prediction 85 / 86 / 92, matched 84 / 77 / 75), the same flat run as the
  README's bank effect (84 / 83 / 87); white measured 33 / 37 / 41 (bank 35 / 36 / 50, matched 35 / 35 / 38). The bank, not the matched template, describes the injections.
  The white values are below the README's theory-only 42 / 43 / 50 m; the thresholds differ (synthetic 2-of-3 at 1/h: 3.5-3.7 against the real-data
  ones), not checked further.
- Statistics are small: R50 is read off a geometric grid with a step of 1.46 and 10-30 injections per point, single points differ by up to 25 % from the prediction.
  A full run (`--runs-per-mag 2 --trials 30`) takes a few hours.
- With each star's own bank (the default since the per-star bank, `synthetic_min_size.py` without `--fixed-bank`), one white B = 9 run
  (`synthetic_smoke_star`, 6 trials x 5 radii, so only indicative): the matched and the bank prediction coincide (2-of-3 R50 39.9 m both, measured 40.0 m)
  and measured / predicted SNR is 0.96-0.99, against 0.84 (matched) with the fixed bank at the same magnitude: the bank loss is gone for a star of known size.
  A full per-star run (white and coloured, B = 7 / 8 / 9) is still to be done.

## Open

- **Per-star search bank: code done 2026-11-16, regeneration pending.** `mimic.star_bank` builds the bank from the star's own size
  (B-magnitude relation) and temperature; `run_comparison`, `run_mimic`, `calibrate_run` use it (`MimicConfig.per_star_size`, CLI
  `--star-mas` fixes it again). Smoke tests (nothing in `outputs/` touched): the bank's expected SNR against the perfect template is
  0.99 / 0.98 / 0.97 at B = 7 / 8 / 9 (fixed bank 0.99 / 0.77 / 0.63; LST-1 0.45 at B = 9); on HD170263 (B 8.62) the injected peak SNR is
  1.3-1.5x higher with the own bank, measured / predicted at R = 150 m is 1.04 (fixed bank 0.69), fitted scale 1.02 (0.75). Still to do:
  regenerate `real_candidates.csv` (`noise_removal_comparison.py`, now with a `star_mas` column; `min_size_calibration.py` refuses a
  list without it unless `--force`), rerun the calibration and the study; the R50 numbers above are still those of the fixed-bank pipeline
  for thresholds and calibration.
- `excess_scaling` ("star" vs "total"), needs a focused / defocused LST-1 run, or more nights with a wide range of star fractions.
- The thresholds are extrapolated beyond the measured tail (about 1 telescope-hour of candidates; 3 nights, only 2 with 3 channels).
- Catalogue of target stars: run on the real catalogue (`examples/studies/catalogue_min_size.py`, top 5 by score without HD 165132, plus dummy
  stars on the relation fitted to them); see README "Catalogue targets". The remaining caveat is the fixed-bank thresholds/calibration.
- The detection width neglects the 100 m telescope offsets, the impact parameter b is not varied inside the calibration, and the
  calibration scatter (1.17) is applied to every b.
- The observing times are per star and per telescope array on-target; no weather, airmass, moon or duty cycle.
- `real_candidates.csv` should be regenerated when nights are added; rerun `min_size_calibration.py` then `min_detectable_size.py`.
