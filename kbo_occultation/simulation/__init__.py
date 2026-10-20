"""
Forward modelling: what an occultation should look like and how bright the
star is for a given telescope.

- engine: Fresnel diffraction light curves (compute_lightcurve, OccultationEngine)
- physics: constants, Planck spectrum, Fresnel/Lommel integrals
- sweep: parameter grids over KBO / star / band with a shared engine
- instruments: telescope throughput, photon counts, noise budget, theoretical SNR
- detectability: spatial -> time conversion, binning, idealised SNR of a dip
- kinematics: shadow velocity of a KBO
"""
