"""Star photon counts: CALSPEC Vega zero point, Johnson B (Vega system), optional colour term."""
import numpy as np
import pytest

from kbo_occultation.analysis.snr import logged_psf_containment
from kbo_occultation.simulation.instruments import Instrument, filter_curve, _trapz
from kbo_occultation.simulation.physics import (H_CGS, T_VEGA_K, VEGA_B_MAGNITUDE, teff_from_spectral_type,
                                                vega_photon_density)
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


@pytest.mark.parametrize("spt, teff", [("B9.5IIIp(HgMn)", 10000), ("B9IV", 10500), ("B5", 15000), ("A0V", 9600)])
def test_teff_from_spectral_type(spt, teff):
    assert teff_from_spectral_type(spt) == pytest.approx(teff)


@pytest.mark.parametrize("spt", [None, "garbage", "", "K2V", "A9"])
def test_teff_unknown_is_none(spt):
    assert teff_from_spectral_type(spt) is None


@pytest.mark.parametrize("tel", ["MAGIC", "LST"])
def test_colour_term(tel):
    inst = Instrument(tel, site="LaPalma")
    assert inst.star_photon_rate(T_VEGA_K) == pytest.approx(inst.star_photon_rate())
    ratios = [inst.star_photon_rate(T) / inst.star_photon_rate() for T in (10000, 12000, 15000)]
    assert ratios == sorted(ratios)
    assert 1.005 <= ratios[-1] / ratios[0] <= 1.07


def test_magic_qe_is_no_filter():
    with pytest.warns(UserWarning, match="MAGIC_QE"):
        inst = Instrument("MAGIC", filter_type="MAGIC_QE", site="LaPalma")
    ref = Instrument("MAGIC", site="LaPalma")
    assert inst.filter_type is None
    assert inst.star_photon_rate() == ref.star_photon_rate()
    assert inst.number_nsb_photons() == ref.number_nsb_photons()


@pytest.mark.parametrize("tel, nsb", [("MAGIC", 0.0977 * 1.103), ("LST", 0.258 * 1.103)])   # hexagonal pixel
def test_nsb_rate_dark_no_filter(tel, nsb):
    assert Instrument(tel, site="LaPalma").number_nsb_photons() == pytest.approx(nsb, rel=0.02)


# ---- noise budget terms -------------------------------------------------------------------------

from kbo_occultation.simulation.instruments import scintillation_variance, pixel_containment  # noqa: E402


def test_hexagon_nsb_is_1p103_times_circle():
    inst = Instrument("MAGIC", site="LaPalma")
    d = inst.pixel_FoV_per_telescope["MAGIC"]
    circle = np.pi / 4 * d ** 2
    assert inst.pixel_solid_angle().value / circle == pytest.approx(1.1027, rel=1e-3)
    inst.pixel_shape_per_telescope["MAGIC"] = "circle"
    assert inst.number_nsb_photons() * 1.1027 == pytest.approx(Instrument("MAGIC", site="LaPalma").number_nsb_photons(), rel=1e-3)


def test_scintillation_regimes():
    D, v = 17.5, 20.0
    long_t = np.array([5.0, 10.0])
    s = scintillation_variance(D, long_t, wind_speed_ms=v)
    assert s[0] / s[1] == pytest.approx(2.0, rel=0.05)                  # t^-1 in variance for t >> D/v
    short = scintillation_variance(D, np.array([1e-4, 1e-5]), wind_speed_ms=v)
    assert short[0] / short[1] == pytest.approx(1.0, rel=0.02)          # saturates for t << D/v
    assert 3e-4 < np.sqrt(short[0]) < 1.2e-3                            # ~7e-4 (MAGIC-sized dish)


def test_pixel_containment_and_psf():
    assert pixel_containment(0.5, 30.0) == pytest.approx(1.0, abs=1e-6)
    lst = Instrument("LST", site="LaPalma")
    assert lst.psf_containment(46.69, "d80") == pytest.approx(0.87, abs=0.02)
    magic = Instrument("MAGIC", site="LaPalma")
    assert 0.6 < magic.psf_containment(9.45, "r39") < 0.8                # sigma = PSF
    assert magic.psf_containment(None, "r39") == 1.0
    with pytest.raises(ValueError):
        magic.psf_containment(9.45, "d39")
    assert logged_psf_containment(magic, 9.45) == magic.psf_containment(9.45, "r39")
    assert logged_psf_containment(lst, 46.69) == lst.psf_containment(46.69, "d80")


def test_psf_scales_star_only_and_pedestal_only_floor():
    tel = Instrument("MAGIC", site="LaPalma")
    a = tel.noise_budget(8.0, 10 * u.ms)
    b = tel.noise_budget(8.0, 10 * u.ms, psf_containment=0.8)
    assert b["star_photons"] == pytest.approx(0.8 * a["star_photons"])
    assert b["nsb_photons"] == pytest.approx(a["nsb_photons"])
    p = tel.noise_budget(8.0, 10 * u.ms, pedestal_pe_per_ns=0.1)
    assert p["sigma_floor"] > a["sigma_floor"]
    for k in ("sigma_shot", "sigma_scint"):
        assert p[k] == pytest.approx(a[k])


def test_excess_noise_is_F_not_F2():
    tel = Instrument("MAGIC", site="LaPalma")
    b = tel.noise_budget(8.0, 10 * u.ms)
    F = tel.excess_noise_factor_per_telescope["MAGIC"]
    assert b["sigma_shot"] == pytest.approx(F * np.sqrt(b["star_photons"] + b["nsb_photons"]) / b["star_photons"])


def test_extinction_k_B_is_imposed():
    tel = Instrument("MAGIC", site="LaPalma")
    tel.set_extinction_k_B(0.23)
    t440 = float(tel.atmospheric_transmission(np.array([440.0]))[0])
    assert -2.5 * np.log10(t440) == pytest.approx(0.23, rel=1e-3)
    tel.set_extinction_k_B(None)
    assert -2.5 * np.log10(float(tel.atmospheric_transmission(np.array([440.0]))[0])) == pytest.approx(0.30, abs=0.02)
