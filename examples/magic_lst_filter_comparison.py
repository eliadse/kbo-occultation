"""
Compare MAGIC vs LST, with vs without the MAGIC_SII filter, to see which
configuration maximizes detectability of a sub-km KBO occultation.

Pipeline per configuration
---------------------------
1. Auto-detect the wavelength range where the instrument (+ filter, if
   any) actually has throughput, and use the instrument's own
   total_transmission / total_transmission_filter as the spectral
   response fed into the diffraction simulation -- so the *shape* of the
   filter (e.g. the SII filter's narrow ~400-440nm passband) is what
   determines which wavelengths contribute to the diffraction pattern,
   not just a flat on/off bandpass.
2. Run a KBO radius x impact-parameter x star-size sweep with that
   response (run_parameter_sweep reuses one OccultationEngine per
   config, so this is cheap).
3. Convert each light curve's spatial axis to time via an assumed shadow
   velocity, and resample to the instrument's real sampling cadence --
   this matters because matched-filter SNR depends on how many
   independent samples the dip spans, which is a property of the real
   DAQ, not of the simulation's spatial grid resolution.
4. Convert the instrument's photon-counting SNR to a per-sample sigma
   (sigma_from_instrument) and compute matched-filter SNR on the
   resampled curve.

READ THIS BEFORE TRUSTING THE FILTERED-VS-UNFILTERED NUMBERS
--------------------------------------------------------------
sigma_from_instrument() -> Instrument.signal_to_noise_ratio() needs the
night-sky-background spectrum (Spectra_NSB_ref.txt /
Spectra_NSB_halfmoon.txt in kbo_occultation/data/iact_reference_values/),
which is NOT currently in the repo, so this script cannot yet compute a
real sigma. Below, PLACEHOLDER_MODE controls what happens:

- PLACEHOLDER_MODE = "flat_nsb": uses Instrument.nsb_per_telescope (a
  flat photons/ns/pixel rate already hardcoded in the class, unused
  elsewhere) instead of the real NSB spectrum. This rate is the SAME
  with or without a filter, so it doesn't capture the filter's main
  benefit: suppressing sky background more than it suppresses your
  narrow-line star signal. That means this mode will systematically
  UNDERSTATE the benefit of using the SII filter (it only sees the
  filter's cost -- fewer star photons -- not its NSB-suppression
  benefit). Treat "no filter wins" results from this mode with real
  suspicion; treat "filter wins anyway" results as trustworthy (they won
  despite the missing NSB benefit, so they'd likely still win, and
  probably by more, with the real spectrum). It's fine for comparing
  MAGIC vs LST *within* the same filter setting, since that comparison
  doesn't depend on the missing file.
- Once Spectra_NSB_ref.txt is added, set PLACEHOLDER_MODE = None to use
  the real Instrument.signal_to_noise_ratio() for everything.
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from kbo_occultation import GridConfig, NumericalConfig, BandpassConfig, run_parameter_sweep
from kbo_occultation.config import standard_sampling, standard_sample_duration
from kbo_occultation.instruments import Instrument
from kbo_occultation.detectability import (
    spatial_to_time,
    resample_to_cadence,
    matched_filter_snr,
    sigma_from_instrument,
)

# ─── Adjust these for your actual target ────────────────────────────
STAR_TEMPERATURE_K = 30000
STAR_MAGNITUDE = 8.0          # apparent magnitude in the instrument's filter/band
SHADOW_VELOCITY_MPS = 25e3     # placeholder; use the real Earth-shadow relative velocity for your event
SITE = "La Palma"              # any value other than 'Paranal' -> La Palma atmospheric transmission
PLACEHOLDER_MODE = None  # set to None once the real NSB spectrum files are added

SAMPLING_DT_S = 20 * standard_sampling * standard_sample_duration * 1e-9  # ~ 0.01s = 10ms

CONFIGS = [
    ("MAGIC", None),
    ("MAGIC", "MAGIC_SII"),
    ("LST", None),
    ("LST", "MAGIC_SII"),
]

RADII_M = [20, 50, 80, 100, 150]
IMPACT_PARAMETERS_M = [0.0, 500.0, 1000.0]
STAR_ANGULAR_RADII_MAS = [0.03]
DISTANCE_AU = 40.0

grid = GridConfig(x_max_m=4000, n_x=800)
numerics = NumericalConfig(n_r_grid=3000)

results = {}
for tel, filt in CONFIGS:
    inst = Instrument(telescope_type=tel, filter_type=filt, site=SITE)
    l_min, l_max, n = inst.response_band()
    band = BandpassConfig(l_min, l_max, n)
    response = inst.total_transmission_filter if filt else inst.total_transmission
    label = f"{tel}{'+SII' if filt else ' (unfiltered)'}"

    band_label = f"{band.lam_min_nm:g}-{band.lam_max_nm:g}nm"
    points = run_parameter_sweep(
        radii_m=RADII_M,
        distances_au=[DISTANCE_AU],
        impact_parameters_m=IMPACT_PARAMETERS_M,
        star_angular_radii_mas=STAR_ANGULAR_RADII_MAS,
        star_temperature_K=STAR_TEMPERATURE_K,
        bands=[band],
        grid=grid,
        numerics=numerics,
        responses={band_label: response},
        keep_curves=True,
    )

    sigma = sigma_from_instrument(inst, STAR_MAGNITUDE, SAMPLING_DT_S * u.s)

    for p in points:
        t_s = spatial_to_time(p.x_m, SHADOW_VELOCITY_MPS)
        t_binned, I_binned = resample_to_cadence(t_s, p.intensity, SAMPLING_DT_S)
        p.snr = matched_filter_snr(I_binned, sigma)
        p.n_samples = len(t_binned)
        p.detectable = p.snr >= 5.0

    results[label] = points
    n_det = sum(p.detectable for p in points)
    print(f"{label:16s} band={band.lam_min_nm:5.1f}-{band.lam_max_nm:5.1f}nm  "
          f"sigma={sigma:.4f}  detectable(SNR>=5): {n_det}/{len(points)}  "
          f"mean SNR={np.mean([p.snr for p in points]):.2f}")

# ─── Plot: SNR vs KBO radius at b=0, one line per configuration ────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

for label, points in results.items():
    subset = sorted(
        (p for p in points if p.impact_parameter_m == 0.0),
        key=lambda p: p.radius_m,
    )
    axes[0].plot([p.radius_m for p in subset], [p.snr for p in subset], "o-", label=label)

# This threshold is a placeholder, I don't know the threshold yet
axes[0].axhline(5.0, color="k", ls="--", lw=1, label="SNR=5 threshold")
axes[0].set_xlabel("KBO radius (m)")
axes[0].set_ylabel("Matched-filter SNR")
axes[0].set_title("Detectability vs KBO radius (b=0)")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)

labels = list(results.keys())
mean_snr = [np.mean([p.snr for p in pts]) for pts in results.values()]
axes[1].bar(labels, mean_snr)
axes[1].set_ylabel("Mean matched-filter SNR across grid")
axes[1].set_title("Overall comparison" + (" (placeholder sigma)" if PLACEHOLDER_MODE else ""))
axes[1].tick_params(axis="x", rotation=20)

plt.tight_layout()
plt.savefig("magic_lst_filter_comparison.png", dpi=200)
plt.show()

best = max(results.items(), key=lambda kv: np.mean([p.snr for p in kv[1]]))
print(f"\nBest configuration by mean SNR: {best[0]}"
      + (" -- re-run with the real NSB spectrum before trusting this if it's an unfiltered config."
         if PLACEHOLDER_MODE else ""))
