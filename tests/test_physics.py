def test_planck_positive():
    import numpy as np
    from kbo_occultation.physics import planck_photon
    
    lam = np.array([500e-9])
    assert planck_photon(lam, 5000) > 0

