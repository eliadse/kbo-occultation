# Noise implementation plan: characterise the diffraction signal, then measure noise that mimics it

Deferred work (noise removal) is in `future-noise-removal.md`.

## Context
- **The noise:** real MAGIC and LST-1 light curves carry 1–4 % RMS of correlated noise at 0.1–1 Hz, plus some 1–10 Hz noise. The likely causes are wind, tracking, vibration and excess scintillation.
- **The goal:** find out how often this noise produces something that looks like an occultation.
- **Two blockers:**
  1. **The signal band is not well established.** A first pass gave f_Nyq of about 13–33 Hz for resolved stars (0.01–0.03 mas, 320–450 nm, 40 AU, 25 km/s). Problems with that pass:
     - 0.02 mas came out lower than 0.03 mas, which is not physical;
     - the numbers depend on a single 99 % power threshold;
     - the grid, the wavelength sampling and the PSD method have not been checked for convergence.
  2. **The search's noise statistics are optimistic.**
     - The local σ is a per-block MAD over ~0.1 s blocks, so it sees only noise above about 10 Hz.
     - The matched filter assumes white noise.
     - Real-data false alarms are counted but have never been characterised.

Point sources are out of scope. Every star is treated as resolved (r★ ≥ 0.01 mas).

## Part 1: Characterise the diffraction signal in frequency (do first)
Files: `kbo_occultation/simulation/detectability.py`, new `kbo_occultation/analysis/bands.py`, new `examples/simulation/signal_band_study.py`.

1. **Audit `spatial_power_spectrum` / `nyquist_frequency`** (`detectability.py:53-95` and the helper they call). Check:
   - **mean handling:** is the PSD taken of I or of 1 − I? A DC term inflates the low-frequency cumulative power;
   - **windowing / zero-padding:** a truncated pattern edge gives spectral leakage;
   - **frequency resolution:** Δf = 1/(2 x_max). At `x_max = 30 km` that is about 0.42 Hz at 25 km/s, and the 99 % point is a single bin with no interpolation, which causes the jumps between bins;
   - **one-sided vs two-sided normalisation.**
   Fix anything wrong. Add an `interpolate=True` option that locates the threshold crossing between bins.
2. **Convergence scan** (each checked against the 0.03 mas baseline):
   - `GridConfig`: x_max of 10, 30 and 60 km; n_x so that dx ≤ 5 m;
   - n_lambda: 11, 41 and 81 across 320–450 nm;
   - `n_star_side`: 16, 32 and 64, for the stellar-disk average (`engine.py:71-116`);
   - `NumericalConfig.n_r_grid`.
   Acceptance: f_Nyq changes by less than 2 % between the two finest settings.
3. **Physical scan** on the converged numerics:
   - r★ = 0.01–0.05 mas in 0.005 mas steps;
   - R = 20–1000 m;
   - impact parameter b = 0–2 km;
   - Δ = 30, 40 and 50 AU;
   - v = 5–30 km/s (quadrature to opposition; use `kinematics`);
   - star temperature 5 000–30 000 K;
   - response: the real `Instrument` MAGIC/LST QE × reflectivity response passed via `response=`, not the dummy trapezoid, over 320–450 nm.
4. **Describe the band with several numbers, not one:**
   - `f50`, `f90`, `f99`: cumulative-power quantiles;
   - `f_low`: the frequency below which 5 % of the power sits, which sets how much a high-pass can remove;
   - the matched-filter-weighted band, i.e. the frequency range holding 90 % of ∫|T(f)|², which equals the filter's sensitivity under white noise;
   - the event chord time H/v from `upper_limits.py:94-108`.
   Check the trends against analytic expectations: f scales with v, and the band narrows as r★/F grows.
5. **Outputs:**
   - `analysis/bands.py::signal_band(star_mas, radii_m, band, distance_au, v, quantiles)` returns these numbers;
   - a figure and a CSV of the scan under `outputs/simulation/signal_band/`;
   - a short README section replacing the 13–33 Hz claim. Update `config.lightcurve_sampling` if the result changes it.
6. **Explain the 0.02 vs 0.03 mas inversion** (a resolution artefact or a real fringe-envelope effect) and record the explanation.

## Part 2: Honest noise normalisation in the search
Files: `kbo_occultation/analysis/matched_filter.py`, `search.py`, `injection_batch.py`.

1. **Empirical filter-output normalisation:**
   - normalise `corr` by a rolling robust σ of `corr` itself (1.4826 × MAD over the detrend window, excluding ±1 template around the peak);
   - implement it as `normalization: Literal["local_block", "global", "empirical"]` (`matched_filter.py:194-206`);
   - default `"empirical"`; keep `"local_block"` for regression.
2. **σ-path test:** compare the `decimate="auto"` and `decimate=1` σ on red-noise surrogates (Part 3.3) and document the gap.
3. **χ² veto:** offer the empirical-σ weighting in `shape_veto_chi2` (`:223-265`). Calibrate `max_chi2_reduced` on off-source real data so that χ²/dof ≈ 1 on pure noise.
4. **Monte Carlo** (`injection_batch.py:367-375`): check on 20 trials whether a per-trial σ of `corr` is needed.

## Part 3: Measure the noise that mimics an occultation
New module: `kbo_occultation/analysis/mimic.py`. CLI subcommand: `kbo mimic` (`cli.py`).

1. **In-band noise per run:**
   - move `segmented_welch` from `examples/studies/*/highpass_filter.py` into `analysis/`;
   - use segments ≥ 2¹⁶ samples when possible, and flag runs where that isn't possible;
   - `inband_rms(lc, band)`: the noise in the Part 1 band;
   - `matched_filter_sigma(lc, template) = √(∫|T|²S / ∫|T|²)`, validated against Part 2's empirical σ (target ≤ 10 % difference);
   - add these columns to `snr_report` (`snr.py`).
2. **Candidate catalogue:**
   - blind search on every logged run (`ObservationLog.select()`) with corrections {raw, dc_detrend, dc_despike} and an SNR threshold of 3;
   - per candidate record: t, channel, SNR, χ²/dof, best R, depth, multiplicity, and the time-slide expectation (`coincidence.py:193-258`);
   - auxiliary features within ±1 s:
     - on-target DC slope and excursion flag (`dc_combine.flag_dc_excursions`);
     - background-pixel DC excursions (a cloud or transparency common mode);
     - starguider deviation and wind (`cc_reports.py`), time-resolved if the reports allow, otherwise as run-level covariates;
   - output to `outputs/mimic/`.
3. **Gaussian surrogates:** phase-randomise each run's spectrum, keeping the gap mask, and draw N = 20 per run. Run the same search. Real false-alarm rate above the surrogate rate points to transient, non-Gaussian mimics such as gusts or glitches.
4. **Negative-tail statistics:** measure the excess over a Gaussian at 3, 5 and 7σ, and its overlap with DC excursions. Nothing is removed; `remove_outliers_lightcurve` deliberately clips only upward spikes.
5. **Deliverable:** FAR(SNR) in events per hour per correction, for single-telescope and coincident events, comparing real data, surrogates and time slides. The table also lists what fraction of the mimics coincide with each auxiliary flag.

## Tests (`tests/unit/`)
- `test_signal_band.py`:
  - f_Nyq scales linearly with v;
  - it is stable to < 2 % under grid refinement;
  - the trend over r★ is monotonic once the artefact is fixed.
- `test_mf_normalization.py`:
  - on white noise, all normalisations agree;
  - on red-noise surrogates, the empirical SNR on pure noise has σ ≈ 1.
- `test_surrogate.py`: the PSD is preserved within Welch errors and the gap mask is kept.

## Verification
1. `pytest tests/unit`.
2. `python examples/simulation/signal_band_study.py`: the convergence plots pass and the band table is written.
3. `kbo mimic --night 2025-12-16`: produces the catalogue and the FAR curves for real data, surrogates and time slides.

## Order
Part 1 → Part 2 → Parts 3.1 and 3.3 → Parts 3.2 and 3.4 → Part 3.5.
