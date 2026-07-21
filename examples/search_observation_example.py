# examples/search_observation_example.py
"""
End-to-end blind coincidence search on one observation -- the entry
point for processing real (on-ecliptic) data, exercised here on a
benchmark-star file.

Loads all three channels (C=MAGIC-1, A=MAGIC-2, B=LST-1), applies a
high-pass correction, runs the per-telescope matched-filter search over
a small grid of template radii, demands >= MIN_NTEL telescopes in
coincidence, and measures the accidental-coincidence background with
time shifts. On a benchmark star (no KBOs expected) any surviving
coincidence is a false positive -- this is how the array pipeline's
false-positive rate is measured.
"""

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.config import BandpassConfig, GridConfig, NumericalConfig, StarConfig
from kbo_occultation.search import search_observation

DATA_DIR = f"{PACKAGE_DATA}/observations"
DATA_FILE = f"{DATA_DIR}/Spectrum_stats_500MSa_Buff_50Ohm_200mV_15min_HD-17316_10002_20002_20251216T223319.bin"
DC_PKL = f"{DATA_DIR}/DCs/2025_12_16/dc_report.pkl"

result = search_observation(
    DATA_FILE,
    # DC reports per telescope (only needed for the DC-based corrections):
    dc_reports={"m1": DC_PKL, "m2": DC_PKL, "lst1": DC_PKL},
    template_radii_m=(100.0, 300.0, 500.0),
    distance_au=43.0,
    star=StarConfig(temperature_K=30000, angular_radius_mas=0.03),
    band=BandpassConfig(300.0, 650.0, 25),
    grid=GridConfig(x_max_m=4000, n_x=800),
    numerics=NumericalConfig(n_r_grid=3000),
    correction="highpass",
    highpass_cutoff_hz=3.0,
    min_ntel=2,
    n_background_shifts=100,
)

print(f"live time: {result.live_time_s:.0f}s, coincidence tolerance: {result.tolerance_s*1e3:.0f}ms")
for ch, res in result.per_telescope.items():
    print(f"  {res.telescope} (ch {ch}): {len(res.candidates)} candidates, sigma={res.sigma:.2e}")

if result.coincidences:
    print(f"\n{len(result.coincidences)} coincidence(s):")
    for ev in result.coincidences:
        members = ", ".join(f"{ch}: SNR={c.snr:.1f} chi2={c.chi2_reduced:.2f}"
                            for ch, c in ev.members.items())
        print(f"  t={ev.time_s:.3f}s n_tel={ev.n_tel} combined SNR={ev.combined_snr:.1f} [{members}]")
else:
    print("\nno coincidences")

if result.accidental is not None:
    print(f"expected accidentals in this observation: "
          f"{result.accidental.expected_in_live_time:.3f} "
          f"(rate {result.accidental.rate_per_s:.2e}/s over {result.accidental.n_shifts} shifts)")
