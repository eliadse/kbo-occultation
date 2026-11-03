# Pending fixes

Open items in the package and the SNR analysis, as of 2026-11-02. Each item lists where the problem
lives and what a fix would do. The completed items of `photon-zeropoint-future-fixes.md` (now deleted) are
recorded in the git history (commits `8161b09`, `7d90d42`).

## Package

1. **LST funnel transmission tuned for a circular pixel** (`simulation/instruments.py`,
   `funnel_transmission_per_telescope['LST'] = 0.892`).
   - The value was tuned to reproduce the reference LST NSB of 0.26516 ph/ns/pixel with a circular pixel.
   - The pixel is now a hexagon (NSB × 1.103), so LST NSB is 0.2844 ph/ns, ~7 % above the value the tuning aimed at.
   - Fix: decide between re-tuning it (~0.81 reproduces 0.2652) and replacing it by a physical funnel transmission.

2. **Electronic pedestal is unmeasured** (`pedestal_pe_per_ns_per_telescope`, 0 for every telescope).
   - It enters only the variance-estimator floor, through ((S + B + P) / S)².
   - The measured NSB is 7–12× the model, and the dark-run NSB in fast units is 0.70–0.89× the dark-run
     mean (DC-to-fast conversion). Both point to a pedestal the DC doesn't see, or to missing star light.
   - Fix: a closed-camera run separates the two. Put the result in the dictionary.

3. **Variance-estimator floor on LST-1** (`variance_estimator_k_per_telescope`).
   - The radiometer equation gives k = 250 MHz / 110 MHz = 2.27 for MAGIC and LST, and it is applied to both.
   - The SNR study saw no floor on LST-1, and LST-1 detr / theo_obs is still ~0.11 at 10 ms (MAGIC ~0.2).
   - Fix: find which is wrong (the floor on LST-1, or something else that makes LST-1 low), then set k for LST.

4. **Other input bands and systems** (`star_photon_rate`).
   - Only Johnson B (Vega system) is supported. A `band` / `system` option (V, or AB bands such as SDSS g or
     Gaia BP) would need the band curve, its system (Vega or AB 3631 Jy) and its Vega magnitude.
   - Not needed while every star has a SIMBAD B magnitude.

5. **Assumed scintillation wind speed** (`scintillation_variance`, `wind_speed_ms=20.0`).
   - It sets the short/long-exposure crossover t ~ D / v and so the short-exposure limit (~7·10⁻⁴ for MAGIC).
   - Fix: replace by a measured value (site turbulence profile or the light curves' scintillation spectrum).

## Notebooks and scripts

6. **`notebooks/snr_review/build_notebook.py` keeps its own copy of the old noise budget.**
   - Scintillation in cm with 1/(2t) (line ~415), the per-bin `sigma_elec` term, the floor on MAGIC only,
     circular pixel, no PSF containment, NSB from the dark-run mean.
   - It is self-contained, so it still runs, but its numbers disagree with `Instrument.noise_budget`.
   - Fix: update it, or make it call `Instrument.noise_budget` and `nsb_table`.

## Analysis

7. **Mean-based vs detrended SNR gap.** Find out why `mean_*` and `detr` differ so much on real data.
8. **High-pass default** (`config.HIGHPASS_CUTOFF_HZ = 3.0`). Re-run
   `examples/studies/noise_removal_comparison.py` with more nights before changing it: the current ranking
   rests on 0.73 telescope-hours.
