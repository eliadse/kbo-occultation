# Verification: star photon zero-point fix

Run after `photon-zeropoint-implementation.md` is done. Report every number asked for, and any failure with its output.

## 1. Unit tests
    pytest tests/unit
All pass, including the new `tests/unit/test_star_photons.py`.

## 2. Zero points and before/after photon counts
In a Python session:
- for MAGIC and LST × (no filter, `MAGIC_SII`), at airmass 1, print:
  - `star_photon_rate()`;
  - the effective Φ0 = rate / (A·1e4·10·∫R dλ);
  - `num_photons_from_mag(1.0, 7.0)`.
- Expected Φ0: ≈ 1144 (no filter) and ≈ 1552 (SII) ph cm⁻² s⁻¹ Å⁻¹.

Expected photo-electrons per ms at B = 7 (no filter):
- MAGIC ≈ 1.08e6 (was 1.91e6);
- LST ≈ 2.75e6 (was 4.85e6).

## 3. Audit notebook agreement
Re-execute `notebooks/snr_theo_study.ipynb`.
- Cell 26: the package curve must coincide with "CALSPEC (Vega colour)" (MAGIC 1.083e6, LST 2.748e6 p.e./ms at B = 7, within 1 %).
- Cell 22: no errors, with legacy values taken from `LEGACY_VEGA`.

## 4. SNR impact table
At 10 ms, B = 6, 8, 10, 12, compute `realistic_signal_to_noise_ratio` for MAGIC and LST × (no filter, SII). It should be close to the "after" column below (before → after, from the audit):

| | B=6 | B=8 | B=10 | B=12 |
|---|---|---|---|---|
| MAGIC no filter | 1079→1065 | 987→910 | 596→425 | 152→90 |
| MAGIC SII | 675→944 | 272→589 | 64→212 | 11→46 |
| LST no filter | 5442→4847 | 3158→2410 | 1105→721 | 241→141 |
| LST SII | 1333→2693 | 443→1093 | 101→341 | 18→72 |

(The "after" values were computed with the flat effective Φ0 of 1144/1552 and 2.512^m. The spectral integral and 10^(0.4 m) may shift them by ~1 %.)

## 5. Real-data SNR report
- Run the SNR report (`kbo-night snr`, or `snr_report(log)` from `kbo_occultation.analysis.snr`).
- Confirm that only `theo_*_B` columns are produced (no `_V`).
- Report the median `ratio_detr_B` per telescope at ~10 ms bins. Compare with the audit values before the fix: MAGIC ≈ 0.12–0.13, LST-1 ≈ 0.023. The ratios should rise by roughly √1.77–1.77×, depending on the S/B regime.

## 6. Filter comparison
- Run `tests/magic_lst_filter_comparison.py` at B = 6, 8, 10, 12 (edit its magnitude constant, or loop).
- Report the SII / no-filter SNR ratio per telescope.
- Expected: no filter still wins, but by ~2–3× (shot noise) instead of ~8×.

## 7. Notebooks and docs
- `python notebooks/snr_review/build_notebook.py` rebuilds without errors, and its star-photon numbers match step 2.
- README renders; the zero-point table in the new subsection matches step 2.

## 8. Memory updates (for the assistant with project memory)
After the numbers above are confirmed, update:
- `theo-snr-audit`: zero point fixed, SII finding;
- `filter-decision-sii-vs-none`: new SII penalty;
- `magic-is-scaled-lst-optics` and `realistic-snr-model`: SNR figures quoted there.
