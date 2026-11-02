# Out of scope for the photon zero-point fix, and future fixes

Items found while auditing the star photon count (2026-10-30) and in `notebooks/snr_theo_study.ipynb`, deliberately
left out of `photon-zeropoint-implementation.md`. Each lists where it lives and what a fix would do.

## Related to the photon count

1. ✅ _Done: `physics.teff_from_spectral_type` + `T_VEGA_K`; `star_photon_rate(temperature_K)` is now a colour term relative to 9600 K (B5: ×1.052, B9: ×1.010); `snr.theoretical_snr_table` passes it and adds `teff_K`. Unknown type = Vega colour_ **Spectral type → temperature colour term.**
   - `Instrument.star_photon_rate(temperature_K=...)` exists, but nothing passes a temperature. `analysis/snr.py` could map `obs.spectral_type` to T_eff (snr_theo_study cell 3 has `SPT_TEFF` and `teff_from_spectral_type`) and pass it to the theory.
   - Effect for the logged B5–B9.5 stars: +0.5 to +6 % star photons.
2. **Other input bands and systems.**
   - Only Johnson B (Vega) is supported. A `band` / `system` option (V, or AB bands such as SDSS g or Gaia BP) would need the band curve, its system (Vega or AB 3631 Jy) and its Vega magnitude.
   - Not needed while every star has a SIMBAD B magnitude.
3. ✅ _Done: `Instrument(filter_type='MAGIC_QE')` warns and is treated as no filter; the file is kept as reference data_ **`filter_type='MAGIC_QE'` double-counts the QE**: the response already contains the QE curve. Remove the file or treat it as a no-op.
4. ✅ _Done: dict removed; MAGIC matched the computed NSB, LST was 2.8 % above; notebook keeps it as `LEGACY_NSB`_ **`nsb_per_telescope` dict is unused** and disagrees with the NSB rates computed by `number_nsb_photons`. Remove it or document it as reference values.
5. ✅ _Done: origin noted in `snr_theo_study.ipynb` and `photon-zeropoint-implementation.md`_ **Origin of the legacy constants** (2020.2 no filter, 356.0 SII, 9.73 gaia_bp) is undocumented. 356 is probably the SDSS u effective wavelength (356 nm) pasted by mistake. If Tarek Hassan's original source turns up, note it in the audit notebook.
6. **Re-run the filter decision.** With the corrected zero points:
   - SII star photons go up 4.4× and no-filter ones down 1.77×;
   - the SII shot-noise SNR penalty shrinks from ~8× to ~2–3×; no filter still wins at B = 6–12 (10 ms).
   - Re-run `tests/magic_lst_filter_comparison.py` and `tests/observing_time_vs_kbo_size.py` at several magnitudes, and revise the conclusion "no filter wins decisively".

## Rest of the theoretical SNR (from the snr_theo_study audit)

7. **Scintillation units** (`Instrument.noise_budget`):
   - Osborn's coefficient is SI (D in metres), but the code uses D in cm (σ ÷ 21.5) and 1/(2t) instead of 1/t (σ ÷ √2), so σ is ~30× too small.
   - The t^(−1/2) law is only valid for t ≫ D/v (~0.5 s). Below that, scintillation saturates at ~6–7·10⁻⁴.
8. **PSF containment** is not modelled.
   - Assumes 100 % of the star light in the pixel.
   - LST-1 (D80 = 46.7 mm) ≈ 87 %.
   - MAGIC: 98–99.7 % if the logged "PSF 39 %" value is a diameter, 62–77 % if it is σ.
9. **Hexagonal pixels**: the NSB uses a circular pixel; a hexagon gives NSB × 1.103.
10. **Variance-estimator floor**:
    - should scale with ((S+B)/S)²;
    - the radiometer equation gives k = 250 MHz / 110 MHz = 2.27, so it should apply to LST-1 too;
    - this conflicts with "LST-1 shows no floor" and is still open.
11. **MAGIC optics** are LST curves × scalars. The 400 nm efficiency is 0.178 vs 0.096 in MAGIC-LST.yml (halo 70 % × PSF 70 %). Distinct MAGIC QE/reflectivity files would make the model and the light-curve shape genuinely MAGIC.
12. **Measured NSB is 7–12× the model**: it may be the electronic pedestal or missing star light. A closed-camera run would separate the two.
13. **Excess-noise factor**: MAGIC-LST.yml calls 1.15 "F⁻¹". Whether it is F or F² changes the shot term by ≤ 15 %.
14. **Electronic noise** is modelled as a 0.7 p.e. readout per bin, which does not fit a variance-based signal.
15. **MODTRAN extinction** (k_B ≈ 0.30) is dustier than the log's k_B = 0.23. Star light is ~6 % low at the zenith.
