"""
Real observations: turn raw files plus the night metadata into LightCurves and
observing conditions.

- log: observation_log.yaml (nights, stars, good-time windows, faulty mirrors)
- prepare: per-night preparation: writes the night's log entry (filenames, data,
  SIMBAD, CC and AMC reports), builds the .npz caches and CC conditions
- simbad: star coordinates / magnitudes / spectral types from SIMBAD
- io: stat-binary reader
- photometry: LightCurve and the .npz pre-processing
- filtering: high-pass and outlier removal
- dc_combine: slow DC reports (loading, combining with, detrending/despiking the fast data)
- cc_reports: MAGIC CC .rep parsing and pointing/moon geometry
"""
