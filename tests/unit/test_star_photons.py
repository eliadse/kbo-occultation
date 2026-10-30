"""Star photon counts: CALSPEC Vega zero point, Johnson B (Vega system), optional colour term."""
import numpy as np
import pytest

from kbo_occultation.simulation.instruments import Instrument, filter_curve, _trapz
from kbo_occultation.simulation.physics import H_CGS, VEGA_B_MAGNITUDE, vega_photon_density
from astropy import units as u


def phi0(tel, filt):
    inst = Instrument(tel, filter_type=filt, site="LaPalma")
    wl = inst.default_wavelengths.astype(float)
    area = inst.mirror_area_per_telescope[tel]
    return inst.star_photon_rate() / (area * 1e4 * 10 * _trapz(inst.response(wl), wl))


@pytest.mark.parametrize("tel", ["MAGIC", "LST"])
def test_zero_point_no_filter(tel):
    assert phi0(tel, None) == pytest.approx(1144, rel=0.03)


@pytest.mark.parametrize("tel", ["MAGIC", "LST"])
def test_zero_point_sii(tel):
    assert phi0(tel, "MAGIC_SII") == pytest.approx(1552, rel=0.03)


def test_martini_b_cross_check():
    wl = np.arange(360.0, 561.0)
    s_b = filter_curve("bessell_B", wl) / wl
    phi = vega_photon_density(wl) * 10 ** (0.4 * VEGA_B_MAGNITUDE)
    assert _trapz(s_b * phi, wl) / _trapz(s_b, wl) == pytest.approx(1392.6, rel=0.04)


def test_ab_vs_vega_b_band():
    wl = np.arange(360.0, 561.0)
    s_b = filter_curve("bessell_B", wl) / wl
    ab = 3631e-23 / (H_CGS * wl * 1e-7) * 1e-8
    b_ab_of_vega = -2.5 * np.log10(_trapz(s_b * vega_photon_density(wl), wl) / _trapz(s_b * ab, wl))
    assert -0.15 <= b_ab_of_vega <= -0.09
    assert b_ab_of_vega - VEGA_B_MAGNITUDE == pytest.approx(-0.15, abs=0.03)


def test_colour_term():
    inst = Instrument("MAGIC", site="LaPalma")
    # no filter: 16 % of the response is below 360 nm, where Vega (Balmer jump) is much fainter than a blackbody
    assert inst.star_photon_rate(9600) == pytest.approx(inst.star_photon_rate(None), rel=0.25)
    assert inst.star_photon_rate(15000) > inst.star_photon_rate(10000)
    # SII lies inside Johnson B: a 9600 K blackbody normalised in B matches Vega there
    sii = Instrument("MAGIC", filter_type="MAGIC_SII", site="LaPalma")
    assert sii.star_photon_rate(9600) == pytest.approx(sii.star_photon_rate(None), rel=0.06)


def test_num_photons_from_mag():
    inst = Instrument("MAGIC", site="LaPalma")
    assert inst.num_photons_from_mag(1.0, 7.0) == pytest.approx(inst.star_photon_rate() * 1e-3 * 10 ** -2.8, rel=1e-12)
    out = inst.num_photons_from_mag(1.0, np.array([6.0, 7.0, 8.0]))
    assert out.shape == (3,) and np.all(np.diff(out) < 0)


def test_noise_budget_temperature():
    inst = Instrument("MAGIC", site="LaPalma")
    default = inst.noise_budget(7.0, 10 * u.ms)
    hot = inst.noise_budget(7.0, 10 * u.ms, temperature_K=12000)
    assert hot["star_photons"] != default["star_photons"]
