"""
Real observations: turn raw files plus the night metadata into LightCurves and
observing conditions.

- log: observation_log.yaml (nights, stars, good-time windows, faulty mirrors)
- prepare: per-night preparation (.bin -> .npz caches, CC conditions, YAML stub)
- io: stat-binary reader
- photometry: LightCurve and the .npz pre-processing
- filtering: high-pass and outlier removal
- dc_combine: slow DC reports (loading, combining with, detrending/despiking the fast data)
- cc_reports: MAGIC CC .rep parsing and pointing/moon geometry
"""
