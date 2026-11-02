# Implementation: star photon count from the CALSPEC Vega spectrum (Johnson B)

Instructions for a coding model working in this repository (`kbo_occultation`). Follow them in order.
Match the surrounding code style (short comments, numpy, astropy). Do not commit.
Verification is a separate file (`photon-zeropoint-verification.md`). Out-of-scope items are in `photon-zeropoint-future-fixes.md`; do not do them.

## Why

`Instrument.num_photons_from_mag` (`kbo_occultation/simulation/instruments.py`) computes

    N = t · A_mirror · 10 · ∫R(λ)dλ · Φ0 · 2.512^(−m)

- R = atmosphere(airmass) × mirror × QE × funnel × telescope transmission × camera window [× optical filter].
- Φ0 is a flat photon density (ph cm⁻² s⁻¹ Å⁻¹ at m = 0) taken from the dict `number_of_photons_from_Vega`:
  - 2020.2 without a filter (key `'QE'`);
  - 356.0 with `MAGIC_SII`.
- m is the star's **Johnson B magnitude** from SIMBAD (`obs.mag_B`), which is a **Vega-system** magnitude.

Folding the CALSPEC Vega spectrum shipped in `kbo_occultation/data/alpha_lyr_stis_011.fits` (scaled to B = 0) with the same R gives:

| case | package Φ0 | correct (Vega) Φ0 | AB 0 mag Φ0 | package / correct |
|---|---|---|---|---|
| no filter, MAGIC & LST | 2020.2 | 1144 | 1289 | 1.77 |
| `MAGIC_SII` | 356.0 | 1552 | 1305 | 0.23 |

Origin of the old constants: 2020.2 (no filter) and 9.73 (gaia_bp) are undocumented; 356.0 is probably the SDSS u effective wavelength (356 nm) pasted in place of a photon density; the ugriz values are AB zero points (Martini's table).

- **Vega system:** m = 0 means "Vega's spectrum". Vega has B = V = 0.03 (B − V = 0), so a B = 0 star is Vega × 10^(0.4·0.03).
- **AB system:** m_AB = −2.5 log10(f_ν / 3631 Jy), i.e. a flat f_ν. In photons, Φ_AB = 5.48·10⁶ / λ[Å] ph cm⁻² s⁻¹ Å⁻¹. In B, B_AB − B_Vega ≈ −0.15 mag.
- Our magnitudes are Vega-system B, so the reference spectrum must be Vega's, not AB.
- The package constants match neither system. The ugriz entries of the dict are AB values from Martini's table, irrelevant here.
- The no-filter and SII responses (effective λ 438 / 420 nm) both lie inside Johnson B, so B is the right input magnitude.

The new model:

    N = t · A · 10 · ∫ R(λ) Φ★(λ) dλ · 10^(−0.4 B)      [λ in nm, ×10 for nm→Å, A in cm²]
    Φ★ = Φ_Vega(λ) · 10^(0.4·0.03)                        (default: Vega-coloured star)
    Φ★ = BB_photon(λ, T) · ∫S′Φ_ref / ∫S′BB_photon        (optional: blackbody normalised to the same B)

Here S′ = S_B(λ)/λ is the Bessell B curve converted from energy- to photon-counting.

Scope decisions (already made):
- B magnitude only (no V, no AB option).
- An optional `temperature_K` colour term defaulting to Vega colours; do NOT map spectral types to temperatures.
- Remove the broken dead code.

## Step 1: Vega reference in `kbo_occultation/simulation/physics.py`

Add under the "Spectral functions" section (`kbo_occultation.PACKAGE_DATA` is defined before any subpackage import, so importing it here is safe):

```python
from functools import lru_cache
from astropy.io import fits
from kbo_occultation import PACKAGE_DATA

H_CGS = 6.62607015e-27     # erg s
C_CGS = 2.99792458e10      # cm s^-1
# Johnson magnitude of Vega: V = 0.03, B - V = 0.00 (Bessell & Murphy 2012), so B = 0.03
VEGA_B_MAGNITUDE = 0.03


@lru_cache(maxsize=1)
def _vega_calspec():
    """CALSPEC alpha_lyr_stis_011: wavelength [A] and photon flux density [ph cm^-2 s^-1 A^-1]."""
    d = fits.getdata(f"{PACKAGE_DATA}/alpha_lyr_stis_011.fits")
    lam_a = np.asarray(d["WAVELENGTH"], dtype=float)
    flam = np.asarray(d["FLUX"], dtype=float)                 # erg s^-1 cm^-2 A^-1
    return lam_a, flam * lam_a * 1e-8 / (H_CGS * C_CGS)       # f_lambda / (hc / lambda)


def vega_photon_density(wl_nm):
    """Vega (CALSPEC) photon flux density above the atmosphere [ph cm^-2 s^-1 A^-1] at wl_nm."""
    lam_a, phot = _vega_calspec()
    return np.interp(np.asarray(wl_nm, dtype=float) * 10.0, lam_a, phot, left=0.0, right=0.0)
```

## Step 2: Johnson B curve data file

Create `kbo_occultation/data/optical_filters/optical_filter_bessell_B.txt`:
- two columns, wavelength [nm] and transmission [%], the same format as the other filter files (the loader divides by 100);
- `#` header lines citing Bessell (1990, PASP 102, 1181) Table 2;
- a header note: "energy-counting response: divide by wavelength for photon-counting integrals (Bessell & Murphy 2012)".

Values (nm, %):

    360 0.0 | 370 3.0 | 380 13.4 | 390 56.7 | 400 92.0 | 410 97.8 | 420 100.0 | 430 97.8 | 440 93.5 | 450 85.3
    460 74.0 | 470 64.0 | 480 53.6 | 490 42.4 | 500 32.5 | 510 23.5 | 520 15.0 | 530 9.5 | 540 4.3 | 550 0.9 | 560 0.0

Check that `astropy.io.ascii.read(..., guess=True)` (as used by `optical_filter_transmission`) still yields `col1`/`col2` with the header comments.

## Step 3: `kbo_occultation/simulation/instruments.py`

1. Imports: `from .physics import vega_photon_density, planck_photon, VEGA_B_MAGNITUDE, nm_m`. Add a module-level `_trapz = getattr(np, "trapezoid", None) or np.trapz` and use it everywhere `trapz` appears in this file.
2. Module-level helper extracted from `optical_filter_transmission`:
   ```python
   def filter_curve(name, wavelength):
       """Transmission [0-1] of data/optical_filters/optical_filter_<name>.txt (in %) on `wavelength` [nm]."""
   ```
   - Same body as the current file-reading part: `ascii.read(..., guess=True)`, `interp1d` inside the file's range, zero outside, `/100`.
   - `optical_filter_transmission` then returns `np.ones_like(wavelength)` for no filter, else `filter_curve(self.filter_type, wavelength)`.
3. New method `response(self, wavelength=None)`:
   - returns `total_transmission_filter(wl)` if `self.filter_type` else `total_transmission(wl)`;
   - this is the star response, atmosphere included;
   - use it in `response_band` too (replace its inline conditional).
4. New method:
   ```python
   def star_photon_rate(self, temperature_K=None):
       """
       Photo-electrons per second from a star of Johnson B = 0 (Vega system), at self.airmass,
       through the full response (atmosphere x telescope x optical filter, if any).

       The star spectrum is CALSPEC Vega scaled to B = 0 (Vega has B = 0.03) unless temperature_K
       is given: then a blackbody of that temperature normalised to the same Bessell B photon count.
       """
   ```
   - `wl = self.default_wavelengths.astype(float)`.
   - `phi_ref = vega_photon_density(wl) * 10 ** (0.4 * VEGA_B_MAGNITUDE)`.
   - If `temperature_K` is not None:
     - `s_b = filter_curve("bessell_B", wl) / wl` (photon-counting);
     - `bb = planck_photon(wl * nm_m, temperature_K)`;
     - `phi = bb * _trapz(s_b * phi_ref, wl) / _trapz(s_b * bb, wl)`.
   - Otherwise `phi = phi_ref`.
   - Return `self.mirror_area_per_telescope[self.telescope_type] * 1e4 * _trapz(self.response(wl) * phi, wl) * 10`, with comments for the unit factors: m²→cm², nm→Å.
5. Replace `num_photons_from_mag`:
   ```python
   def num_photons_from_mag(self, obs_times, mag, temperature_K=None):
       """
       Star photo-electrons in obs_times [ms] for a star of Johnson B magnitude `mag` (Vega system),
       whatever the optical filter. See star_photon_rate.
       """
       return obs_times * 1e-3 * self.star_photon_rate(temperature_K) * 10 ** (-0.4 * np.asarray(mag))
   ```
   `10**(-0.4 m)` replaces `2.512**(-m)`. The difference is 0.0004 mag per mag (0.3 % at B = 8), which is intended.
6. Add a `temperature_K=None` keyword to `signal_to_noise_ratio` and `noise_budget`, passed to `num_photons_from_mag`.
   - `realistic_signal_to_noise_ratio` already forwards `**kwargs`; mention `temperature_K` in its docstring.
   - Document it in `noise_budget`'s Parameters: "Optional blackbody colour of the star (default: Vega colours, B − V = 0)."
7. Docstrings: wherever `magnitude` is described ("Apparent magnitude in the current filter (unfiltered if None)"), change it to "Johnson B magnitude of the star (Vega system), whatever `filter_type` is". In `noise_budget`'s Returns, say `star_photons` are photo-electrons (the response includes the QE).
8. Delete:
   - the `number_of_photons_from_Vega` dict, its Martini-URL comment, and its line in the class docstring ("Loaded data");
   - the `# TODO Move the Vega stuff to physics ...` comment;
   - `star_flux` and `star_flux_folded_with_transmission` (dead; the latter calls an unbound `star_flux`).
9. Add one comment line above `qes_per_telescope` or near `optical_filter_transmission`: "`filter_type='MAGIC_QE'` multiplies the QE curve a second time (the response already includes the QE)."

## Step 4: callers

- `kbo_occultation/analysis/snr.py`:
  - set `BANDS = ("B",)`;
  - in the module docstring (around lines 32–35), say the theory uses the B magnitude with the CALSPEC Vega zero point, and remove any mention of V theory columns.
  - Nothing else changes: `theoretical_snr_table` loops over `bands`, and `plot_snr_vs_magnitude` already defaults to `("B",)`.
  - `mag_V` stays in the observation log and `Observation.magnitude("V")` keeps working.
- `kbo_occultation/simulation/detectability.py` `sigma_from_instrument`:
  - docstring `magnitude`: "Johnson B magnitude (Vega system)";
  - mention `temperature_K` among the `**noise_kw` forwarded to `noise_budget`;
  - on the `realistic=False` branch, pass `temperature_K=noise_kw.get("temperature_K")` to `signal_to_noise_ratio`.
- `kbo_occultation/analysis/noise_removal.py` (~line 274): already passes `obs.mag_B`; no change.
- Scripts `tests/snr_vs_magnitude.py`, `tests/magic_lst_filter_comparison.py`, `tests/observing_time_vs_kbo_size.py`, `examples/simulation/plot_lightcurve_comparison.py`:
  - no code change;
  - fix comments or docstrings that say the magnitude is "in the filter band" or "same input mag for both bands" to say Johnson B for every filter.
- `notebooks/snr_review/build_notebook.py`:
  - remove `VEGA_PHOTONS_QE = 2020.2` and make `star_photons(kind, magnitude, bin_ms, airmass)` use the package (set the airmass on the `Instrument` around `num_photons_from_mag`, as `noise_budget(extinction=True)` does);
  - rewrite the markdown around lines 332–336: "S = photo-electrons of a star of Johnson B magnitude m: CALSPEC Vega spectrum × response, scaled by 10^(−0.4(m − 0.03))";
  - rebuild the notebook with the script.
- `notebooks/snr_theo_study.ipynb` (untracked audit notebook):
  - cells 20, 22, 26 and 47 read `INSTR[...].number_of_photons_from_Vega`. Add `LEGACY_VEGA = {'u': 1539.3, 'g': 1134.6, 'r': 875.4, 'i': 714.5, 'z': 602.2, 'MAGIC_SII': 356.0, 'QE': 2020.2, 'MAGIC_QE': 2020.2, 'gaia_bp': 9.73}` to cell 3 and use it there.
  - Where cell 26 compares "package" photons, label the old curve "legacy 2020.2". The package curve should now coincide with "CALSPEC (Vega colour)".
  - In the verdict cell 48, mark item 1 as fixed and add the SII finding (356 → 1552, 4.4× too low) and the AB/Vega numbers from the table above.
  - Edit with NotebookEdit or a small JSON script, then re-execute.

## Step 5: tests (new `tests/unit/test_star_photons.py`)

Use `Instrument(tel, filter_type=f, site="LaPalma")`, with `wl = inst.default_wavelengths.astype(float)` and `phi0 = rate / (A * 1e4 * 10 * trapz(inst.response(wl), wl))`.

1. No filter, MAGIC and LST: `phi0 == approx(1144, rel=0.03)`.
2. `MAGIC_SII`, MAGIC and LST: `phi0 == approx(1552, rel=0.03)`.
3. Martini cross-check: the photon-counting Bessell B weighted mean of `vega_photon_density × 10^(0.4·0.03)` is `approx(1392.6, rel=0.04)`.
4. AB documentation test:
   - define, in the test only, `ab = 3631e-23 / (H_CGS * wl_nm * 1e-7) * 1e-8`;
   - Vega's B_AB = −2.5 log10(∫S′·vega / ∫S′·ab) must be in [−0.15, −0.09];
   - (B_AB − B_Vega) = B_AB − 0.03 ≈ −0.15 ± 0.03.
5. Colour term:
   - `star_photon_rate(9600)` is within 2 % of `star_photon_rate(None)` for MAGIC, no filter;
   - `star_photon_rate(15000) > star_photon_rate(10000)`.
6. `num_photons_from_mag(1.0, 7.0)` equals `star_photon_rate() * 1e-3 * 10**(-2.8)` (rel 1e-12). It accepts an array of magnitudes.
7. `noise_budget(7.0, 10*u.ms, temperature_K=12000)` runs and its `star_photons` differs from the default.

Run `pytest tests/unit`. Existing tests use ratios (test_snr_measurement, test_noise_removal) or only `theo_obs_B` (test_snr_report) and should pass unchanged. If one asserted an absolute photon count, update the expected value and say so in your report.

## Step 6: documentation

- `README.md`, section "Configuration reference": add a subsection **Star photons and magnitude system** (~15 lines) covering:
  - the formula above;
  - the CALSPEC source (`alpha_lyr_stis_011.fits`);
  - "magnitudes are Johnson B, Vega system (Vega B = 0.03)";
  - a two-sentence explanation of AB and why it is not used (B_AB − B_Vega ≈ −0.15 mag; different reference shape: AB/Vega photons = 1.13 with no filter, 0.84 with SII);
  - why B fits both the no-filter and SII responses;
  - the optional `temperature_K` colour term;
  - a small table of the B = 0 photo-electron rate per second and effective Φ0 for MAGIC/LST × no filter/SII, filled from `star_photon_rate()` at airmass 1.
- Update the "Bundled data" paragraph to mention `optical_filter_bessell_B.txt` and that `alpha_lyr_stis_011.fits` is now used.

## Report back

List the files changed, the before/after star photo-electrons per ms at B = 7 for MAGIC/LST × no filter/SII, the test results, and anything you could not do.
