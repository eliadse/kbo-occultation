"""
Functions and reference files created and compiled by Tarek Hassan
"""
import warnings

import numpy as np
from astropy.io import ascii
from astropy import units as u
from scipy import interpolate

from kbo_occultation import PACKAGE_DATA
from .physics import vega_photon_density, planck_photon, VEGA_B_MAGNITUDE, T_VEGA_K, nm_m

_trapz = getattr(np, "trapezoid", None) or np.trapz  # np.trapz removed in NumPy 2.0

# one fast-photometry frame: 2**18 samples x 2 ns
FRAME_SAMPLES_DURATION_S = 262144 * 2e-9
reference_files_path = "{}/iact_reference_values".format(PACKAGE_DATA)

H_SCALE_M = 8000.0                       # atmospheric scale height in the scintillation formulas
HEX_OVER_CIRCLE = (np.sqrt(3) / 2) / (np.pi / 4)   # hexagon (flat-to-flat d) / circle (diameter d) area = 1.103
R80_OVER_SIGMA = np.sqrt(2 * np.log(5))  # 2D Gaussian: r80 = 1.794 sigma; r39 = sigma
# Variance-estimator k from the radiometer equation: Nyquist of the 2 ns samples over the
# signal bandwidth b_v of MAGIC-LST.yml (250 MHz / 110 MHz).
RADIOMETER_K = 250e6 / 110e6


def scintillation_variance(D_m, t_s, airmass=1.0, altitude_m=2200.0, coeff=1.5, wind_speed_ms=20.0):
    """
    Fractional scintillation variance of a star seen through an aperture of diameter D_m [m] in an
    exposure of t_s [s] (Young 1967, in the form of Osborn et al. 2015 eq. 7, SI units).

    Long exposure (t >> D/v):  sigma^2 = 1e-5 C^2 D^(-4/3) X^3 exp(-2h/H) / t
    Short exposure (t << D/v): sigma^2 = 17.34 D^(-7/3) X^3 W   (Kenyon et al. 2006), independent of t,
    with W = int h^2 Cn^2 dh fixed so that both forms agree at the wind speed v.
    The two are joined as 1/sigma^2 = 1/sigma_short^2 + 1/sigma_long^2.
    """
    t_s = np.asarray(t_s, dtype=float)
    atm = airmass ** 3 * np.exp(-2.0 * altitude_m / H_SCALE_M)
    long_var = 1e-5 * coeff ** 2 * D_m ** (-4.0 / 3.0) * atm / t_s
    W = 1e-5 * coeff ** 2 * np.exp(-2.0 * altitude_m / H_SCALE_M) / 10.66 * wind_speed_ms
    short_var = 17.34 * D_m ** (-7.0 / 3.0) * airmass ** 3 * W
    return long_var * short_var / (long_var + short_var)


def pixel_containment(sigma_mm, flat_mm, shape='hexagon', n=601):
    """Fraction of a centred circular Gaussian (sigma_mm) inside a pixel of flat-to-flat size flat_mm."""
    a = flat_mm / 2.0
    half = a * 2.0 / np.sqrt(3) if shape == 'hexagon' else a
    x = np.linspace(-half, half, n)
    X, Y = np.meshgrid(x, x)
    inside = np.ones_like(X, dtype=bool)
    angles = [0, 60, 120] if shape == 'hexagon' else [0, 90]
    for ang in np.deg2rad(angles):
        inside &= np.abs(X * np.cos(ang) + Y * np.sin(ang)) <= a
    g = np.exp(-(X ** 2 + Y ** 2) / (2 * sigma_mm ** 2)) / (2 * np.pi * sigma_mm ** 2)
    return float(g[inside].sum() * (x[1] - x[0]) ** 2)

def load_response_file(filename):
    """
    Load wavelength-dependent response from a text 
    file stored inside kbo_occultation/data/

    Expected format:
    wavelength_nm   transmission_percent
    """

    data = np.loadtxt("{}/{}".format(PACKAGE_DATA, filename))

    lam_nm = data[:, 0]
    response = data[:, 1] / 100.0  # convert % → [0,1]

    return lam_nm, response

def build_response_function(lam_nm_data, response_data):
    """
    Returns an interpolation function R(lambda_nm).
    """

    def response(lam_nm):
        return np.interp(
            lam_nm,
            lam_nm_data,
            response_data,
            left=0.0,
            right=0.0
        )

    return response

def combine_responses(*responses):
    """
    Multiply multiple response functions together.
    """

    def combined(lam_nm):
        R = 1.0
        for r in responses:
            R *= r(lam_nm)
        return R

    return combined

def filter_curve(name, wavelength):
    """Transmission [0-1] of data/optical_filters/optical_filter_<name>.txt (in %) on `wavelength` [nm]."""
    # The guess=True addition let's us use any filter file
    data = ascii.read("{}/optical_filters/optical_filter_{}.txt".format(PACKAGE_DATA, name),
                    guess=True)
    f = interpolate.interp1d(data['col1'], data['col2'])
    result = np.zeros_like(wavelength, dtype=float)
    mask = np.logical_and(wavelength >= np.min(data['col1']), wavelength <= np.max(data['col1']))
    result[mask] = f(wavelength[mask])
    return result/100.

#TODO Add a function that lists all available filters and/or transmission files
class Instrument():
    """

    A class to compute the transmission of the telescopes and their signal to noise ratio

    """
    def __init__(self, telescope_type='MAGIC', filter_type=None, site='Paranal'):
        self.telescope_type = telescope_type
        self.filter_type = self._normalise_filter(filter_type)
        self.site = site

        """
        Parameters
        ----------
        wavelength: float
            Wavelength in nm to calculate the transmission

        telescope_type: str
            Name of telescope. Set to MAGIC.
            Options are LST, MST-N, pSCT, SST-1M, 
            VERITAS, VERITAS++, MAGIC.
        
        filter_type: str
            Type of optical filter. Default is None
            Set to standard filters like g, u, r...
            Add "super_" to select hiperCAM filters
            or "MAGIC_" for the MAGIC filters like SII

        site: str
            Telescope site. Set to Paranal
            Selecting any other selects La Palma stats
        

        Loaded data
        -----------
        default_wavelengths: An array from 100 to 1000 nm

        The night-sky background rate [photons/ns/pixel] is not stored: it is
        computed from the NSB spectrum by number_nsb_photons.

        pixel_FoV_per_telescope: Telescope pixel FoV (diameter)[deg]
        
        mirror_reflectivity_per_telescope: Telescope shadowing?
        
        qes_per_telescope: Quantum efficiency of the telescopes
        
        transmission_per_telescope: Camera elements transmission (plexiglass, etc.)
        
        funnel_transmission_per_telescope: Funnel efficiency at 0 offset
        
        camera_transmission_per_telescope: Camera transmission of the telescopes
        
        mirror_area_per_telescope: Geometric mirror area BEFORE SHADOWING [m^2]

        """
        self.default_wavelengths = np.arange(100, 1000, 1)

        self.pixel_FoV_per_telescope = {
            'MAGIC': 0.1,
            'LST': 0.102,
            'MST-N': 0.179,
            'pSCT': 0.069,
            'SST-1M': 0.241,
            'VERITAS': 0.15,
            'VERITAS++': 0.07
            }
        # SST-GCT           0.161 deg
        # SST-ASTRI         0.188 deg

        # Pixel shape: the FoV above is the flat-to-flat size. The PMT cameras have hexagonal
        # pixels (area x 1.103 relative to a circle of that diameter), the pSCT SiPMs are square.
        self.pixel_shape_per_telescope = {
            'MAGIC': 'hexagon', 'LST': 'hexagon', 'MST-N': 'hexagon', 'pSCT': 'square',
            'SST-1M': 'hexagon', 'VERITAS': 'hexagon', 'VERITAS++': 'hexagon'
            }

        # Focal length [mm], to convert the pixel FoV and the logged PSF to the camera plane.
        self.focal_length_mm_per_telescope = {'MAGIC': 17000.0, 'LST': 28000.0, 'MST-N': 16000.0}

        # Telescope optical PSF requirement FoV:
        # LST               PSF θ_80 of <0.11 degrees, up to 1.2 degrees
        # MST               PSF θ 80 of <0.18 degrees, up to 2.8 degrees
        # MST-SCT           PSF θ 80 of <0.18 degrees, up to 2.8 degrees
        # SST-1M            PSF θ 80 of <0.25 degrees, up to 3.2 degrees
        # SST-GCT           PSF θ 80 of <0.25 degrees, up to 3.2 degrees
        # SST-ASTRI         PSF θ 80 of <0.25 degrees, up to 3.2 degrees

        # Telescope shadowing
        self.mirror_reflectivity_per_telescope = {
            'MAGIC': 'LST_ref_200_1100_190211a.dat',
            'LST': 'LST_ref_200_1100_190211a.dat',
            'MST-N': 'MST_ref_AlSiO2HfO2.dat',
            'pSCT': 'pSCT_Reflectance_SC-MST_Prod3.dat',
            'SST-1M': 'SST-1M_ref_SST1M_DOJLO_ave_extrapol.dat',
            'VERITAS': 'LST_ref_200_1100_190211a.dat',
            'VERITAS++': 'LST_ref_200_1100_190211a.dat'
            }

        self.qes_per_telescope = {
            'MAGIC': 'LST_qe_R11920-100-02.dat',
            'LST': 'LST_qe_R11920-100-02.dat',
            'MST-N': 'MST-N_qe_R12992-100-05.dat',
            'pSCT': 'pSCT_PDE_SiPM_SC-MST_Prod3x.dat',
            'SST-1M': 'SST-1m_PDE_V_4.4V_LVR5_ext.txt',
            'VERITAS': 'LST_qe_R11920-100-02.dat',
            'VERITAS++': 'LST_qe_R11920-100-02.dat'
            }

        # Transmission per telescope ('telescope_transmission' in simtel_array)
        self.transmission_per_telescope = {
            'MAGIC': 0.8928,
            'LST': 0.9606,
            'MST-N': 0.8928,
            'pSCT': 0.81,
            'SST-1M': 0.82,
            'VERITAS': 0.8928,
            'VERITAS++': 0.8928
            }

        self.funnel_transmission_per_telescope = {
            'MAGIC': 0.85,
            'LST': 0.892,        # tweaked to match the reference LST NSB rate 0.26516 ph/ns/pixel with a circular pixel; the pixel is now a hexagon (NSB x 1.103), value left unchanged
            'MST-N': 0.872,
            'pSCT': 1.0,
            'SST-1M': 0.9,
            'VERITAS': 0.85,
            'VERITAS++': 0.85
            }

        self.camera_transmission_per_telescope = {
            'MAGIC': 'MST-N_No7-10_ave.dat',
            'LST': 'LST_transmission_lst_window_No7-10_ave.dat',
            'MST-N': 'MST-N_No7-10_ave.dat',
            'pSCT': 'SST-1M_filter_SST1M_fusedsilica.dat',
            'SST-1M': 'SST-1M_filter_SST1M_fusedsilica.dat',
            'VERITAS': 'MST-N_No7-10_ave.dat',
            'VERITAS++': 'MST-N_No7-10_ave.dat'
            }

        # If no source is specified, the number is coming from here: https://www.mpi-hd.mpg.de/hfm/CTA/MC/Prod4/Config/Efficiencies/
        self.mirror_area_per_telescope = {
            'MAGIC': 240,
            'LST': 380.271, # on-axis, before shadowing
        #    'LST': 391.0, #   Total mirror surface area: 391.0 m**2 from https://forge.in2p3.fr/projects/cta/repository/entry/ASWG/Simulations/MCModelDescription/trunk/datFiles/mirror_CTA-LST-flen_grouped.dat
            'MST-N': 99.03300, # on-axis, before shadowing
        #    'MST-N': 107.2, # Total mirror surface area: 107.2 m**2 from https://forge.in2p3.fr/projects/cta/repository/entry/ASWG/Simulations/MCModelDescription/trunk/datFiles/mirror_CTA-100_1.20-86-0.04.dat
            'pSCT': 50.281, # on-axis, before shadowing
        #    'pSCT': 58.234738, # Coming from primary diameter = 966.38 [cm] and hole 438.66 [cm]
            'ASTRI': 7.931, # on-axis, before shadowing
            'SST-1M': 9.12, # on-axis, before shadowing
        #    'SST-1M': 9.5, # Total mirror surface area: 9.5 m**2 from mirror_DC4m_78cm_s2cm_f5.600_nocent.dat
            'VERITAS': 100,
            'VERITAS++': 100
            }
        
        # Photodetector excess-noise factor F (variance is inflated by F^2
        # relative to pure Poisson shot noise). Values from MAGIC-LST.yml, where
        # the entry is labelled "F^-1"; it is taken to be F itself (shot variance x F^2 = 1.32).
        self.excess_noise_factor_per_telescope = {
            'MAGIC': 1.15,
            'LST': 1.15,
            'MST-N': 1.15,
            'pSCT': 1.15,
            'SST-1M': 1.15,
            'VERITAS': 1.15,
            'VERITAS++': 1.15
            }

        # Electronic pedestal: photoelectron-equivalent background per ns that is present with no
        # light (the variance of the electronics). It adds to the variance the estimator measures
        # but not to the star signal, so it only enters the floor, through ((S + B + P) / S)^2.
        # Unmeasured (0): a closed-camera run would give it (see pending-fixes.md, item 2).
        self.pedestal_pe_per_ns_per_telescope = {t: 0.0 for t in self.excess_noise_factor_per_telescope}

        # Variance-estimator floor: the fast signal is a variance (std^2) of N
        # digitizer samples per frame, whose relative error is sqrt(2k/N) with
        # k the number of effectively independent samples. Radiometer equation:
        # k = Nyquist / signal bandwidth = 250 MHz / 110 MHz (b_v in MAGIC-LST.yml) = 2.27,
        # next to the k ~ 2 measured on MAGIC (N = 2^18; see the SNR study). It applies to
        # LST-1 too, although the SNR study saw no floor there (open question).
        self.variance_estimator_k_per_telescope = {t: 0.0 for t in self.excess_noise_factor_per_telescope}
        self.variance_estimator_k_per_telescope['MAGIC'] = RADIOMETER_K
        self.variance_estimator_k_per_telescope['LST'] = RADIOMETER_K
        self.variance_estimator_n_samples = 2 ** 18

        # Observatory altitude (m) for the atmospheric scintillation term.
        # Paranal ~2635 m, ORM/La Palma ~2200 m.
        self.site_altitude_m = 2635.0 if self.site == 'Paranal' else 2200.0

        # Airmass at which the atmospheric extinction is evaluated (1 = zenith).
        # See set_airmass / noise_budget(extinction=True).
        self.airmass = 1.0

        # B-band extinction [mag/airmass] to impose on the MODTRAN optical depth, or None to
        # keep MODTRAN as it is (k_B = 0.30 at La Palma). See set_extinction_k_B.
        self.extinction_k_B = None

        self.transmission = self.total_telescope_transmission()

    def print_parameters(self):
        print(f"Telescope type: {self.telescope_type}")
        print(f"Optical filter: {self.filter_type}")
        print(f"Telescope site: {self.site}")

    def set_airmass(self, airmass):
        """
        Evaluate the atmospheric extinction (and hence every star photon count)
        at this airmass instead of at the zenith: transmission = exp(-tau * X).
        """
        self.airmass = float(airmass)

    def set_extinction_k_B(self, k_B):
        """
        Scale the MODTRAN optical depth so that the B-band (440 nm) extinction is k_B [mag/airmass]
        (e.g. the observation log's 0.23, against MODTRAN's 0.30), keeping the MODTRAN wavelength
        dependence. None restores MODTRAN.
        """
        self.extinction_k_B = None if k_B is None else float(k_B)

    @staticmethod
    def _normalise_filter(filter_type):
        """'MAGIC_QE' is the QE curve, which the response already contains: treat it as no filter."""
        if filter_type == 'MAGIC_QE':
            warnings.warn("filter_type='MAGIC_QE': the telescope response already includes the QE curve, "
                          "treating it as no filter", stacklevel=3)
            return None
        return filter_type

    def set_filter(self, filter_type):
        """
        Set or change the optical filter of the telescope.
        """
        self.filter_type = self._normalise_filter(filter_type)

    def mirror_reflectivity(self, wavelength=None):
        if wavelength is None:
            wavelength = self.default_wavelengths
        # Interpolate mirror reflectivity for a given IACT:
        data = ascii.read("{}/{}".format(reference_files_path, self.mirror_reflectivity_per_telescope[self.telescope_type]),
                        data_start=1, guess=False, delimiter=',', format='no_header')
        f = interpolate.interp1d(data['col1'], data['col2'])
        result = np.zeros_like(wavelength, dtype=float)
        mask = np.logical_and(wavelength >= np.min(data['col1']), wavelength <= np.max(data['col1']))
        result[mask] = f(wavelength[mask])
        if self.telescope_type == 'VERITAS' or self.telescope_type == 'VERITAS++':
            result = result * 0.6
        elif self.telescope_type == 'MAGIC':
            result = result * 0.82
        return result

    def quantum_efficiency(self, wavelength=None):
        if wavelength is None:
            wavelength = self.default_wavelengths
        # Interpolate quantum efficiency for a given telescope:
        data = ascii.read("{}/{}".format(reference_files_path, self.qes_per_telescope[self.telescope_type]),
                        data_start=1, guess=False, delimiter=',', format='no_header')
        f = interpolate.interp1d(data['col1'], data['col2'])
        result = np.zeros_like(wavelength, dtype=float)
        # The range of wavelengths is wide enough to cover all possible filters
        # so for each one we get their lower and upper bounds and do the calculations
        # within those bounds to be faster.
        mask = np.logical_and(wavelength >= np.min(data['col1']), wavelength <= np.max(data['col1']))
        result[mask] = f(wavelength[mask])
        if self.telescope_type == 'VERITAS' or self.telescope_type == 'VERITAS++':
            result = result * 0.85
        elif self.telescope_type == 'MAGIC':
            result = result * 0.86
        return result

    def camera_transmission(self, wavelength=None):
        if wavelength is None:
            wavelength = self.default_wavelengths
        # Interpolate camera transmission for a given IACT:
        data = ascii.read("{}/{}".format(reference_files_path, self.camera_transmission_per_telescope[self.telescope_type]),
                        data_start=1, guess=False, delimiter=',', format='no_header')
        f = interpolate.interp1d(data['col1'], data['col2'])
        result = np.zeros_like(wavelength, dtype=float)
        mask = np.logical_and(wavelength >= np.min(data['col1']), wavelength <= np.max(data['col1']))
        if self.telescope_type == 'MST-N' or self.telescope_type == 'VERITAS' or self.telescope_type == 'VERITAS++' or self.telescope_type == 'MAGIC':
            result[mask] = f(wavelength[mask])/100.
        else:
            result[mask] = f(wavelength[mask])
        return result

    def total_telescope_transmission(self, wavelength=None):
        """
        Compute the total transmission of the telescope in a wavelength range.
        Takes into account mirror reflectivity, quantum efficiency, funnel
        transmission, camera elements transmission and camera transmission.

        Parameters
        ----------
        wavelength : 
            Array of wavelengths in nm. 
            If set to None, wavelength = np.arange(100, 1000, 1)

        """
        if wavelength is None:
            wavelength = self.default_wavelengths
        resulting_transmission = self.mirror_reflectivity(wavelength) * \
                                self.quantum_efficiency(wavelength) * \
                                self.transmission_per_telescope[self.telescope_type] * \
                                self.funnel_transmission_per_telescope[self.telescope_type] * \
                                self.camera_transmission(wavelength)
        return resulting_transmission

    def atmospheric_transmission(self, wavelength=None):
        """
        Compute the transmission of the atmosphere at a given site in a wavelength range.
        
        Parameters
        ----------
        wavelength : 
            Array of wavelengths in nm. 
            If set to None, wavelength = np.arange(100, 1000, 1)
        
        """
        if wavelength is None:
            wavelength = self.default_wavelengths
        # Interpolate atmospheric transmission for a given CTA site:
        # This function extracts the transmission above 100km altitude
        filename = "atm_trans_2150_1_10_0_0_2150.dat" if self.site == "Paranal" else "atm_trans_2147_1_3_0_0_0.dat"

        col_names = ["wavelength(nm)", "2.197", "2.247", "2.347", "2.447",
                    "2.647", "2.847", "3.147", "3.647", "4.147", "4.500", "5.000",
                    "5.500", "6.000", "7.000", "8.000", "9.000", "10.000", "11.000",
                    "12.000", "13.000", "14.000", "15.000", "16.000", "18.000", "20.000",
                    "22.000", "24.000", "26.000", "28.000", "30.000", "32.500", "35.000",
                    "37.500", "40.000", "45.000", "50.000", "60.000", "70.000", "80.000",
                    "100.000"]

        data = ascii.read("{}/{}".format(reference_files_path, filename), format='no_header',
                        names=col_names, data_start=0, guess=False, delimiter=' ')
        f = interpolate.interp1d(data['wavelength(nm)'], data["100.000"])
        result = np.zeros_like(wavelength, dtype=float)
        mask = np.logical_and(wavelength >= np.min(data['wavelength(nm)']), wavelength <= np.max(data['wavelength(nm)']))
        scale = 1.0
        if self.extinction_k_B is not None:
            tau_b = float(f(440.0))
            scale = self.extinction_k_B / (2.5 * np.log10(np.e) * tau_b)
        result[mask] = np.exp(-self.airmass * scale * f(wavelength[mask]))
        return result

    def total_transmission(self, wavelength=None):
        """
        Combined transmission of the atmosphere and telescope without filters.
        
        Parameters
        ----------
        wavelength : 
            Array of wavelengths in nm. 
            If set to None, wavelength = np.arange(100, 1000, 1)
        
        """
        if wavelength is None:
            wavelength = self.default_wavelengths
        tot_transmission = self.atmospheric_transmission(wavelength) * \
                        self.total_telescope_transmission(wavelength)
        return tot_transmission

    def total_transmission_filter(self, wavelength=None):
        """
        Combined transmission of the atmosphere, the telescope 
        and the chosen filter
                
        Parameters
        ----------
        wavelength : 
            Array of wavelengths in nm. 
            If set to None, wavelength = np.arange(100, 1000, 1)
        """
        if wavelength is None:
            wavelength = self.default_wavelengths
        tot_transmission = self.atmospheric_transmission(wavelength) * \
                        self.total_telescope_transmission(wavelength) * \
                        self.optical_filter_transmission(wavelength)
        return tot_transmission
        
    def optical_filter_transmission(self, wavelength=None):
        """
        Compute the transmission of the chosen optical filter in the wavelength range

        Parameters
        ----------
        wavelength : 
            Array of wavelengths in nm. 
            If set to None, wavelength = np.arange(100, 1000, 1)
        
        """
        if wavelength is None:
            wavelength = self.default_wavelengths
        # If no filter given, return just ones:
        if self.filter_type is None:
            return np.ones_like(wavelength)
        return filter_curve(self.filter_type, wavelength)

    def response(self, wavelength=None):
        """Star response: atmosphere x telescope x optical filter (if any)."""
        if self.filter_type:
            return self.total_transmission_filter(wavelength)
        return self.total_transmission(wavelength)

    def response_band(self, threshold=0.01, margin_nm=10.0, n_lambda=60):
        """
        Auto-detect the useful wavelength range for an instrument+filter
        combo (where its total throughput is above `threshold` of its peak),
        instead of guessing a bandpass by hand.
        """
        wl = self.default_wavelengths.astype(float)
        resp = self.response(wl)
        peak = resp.max()
        above = wl[resp > threshold * peak]
        lam_min = max(wl.min(), above.min() - margin_nm)
        lam_max = min(wl.max(), above.max() + margin_nm)
        return float(lam_min), float(lam_max), n_lambda

    def star_photon_rate(self, temperature_K=None):
        """
        Photo-electrons per second from a star of Johnson B = 0 (Vega system), at self.airmass,
        through the full response (atmosphere x telescope x optical filter, if any).

        The star spectrum is CALSPEC Vega scaled to B = 0 (Vega has B = 0.03). If temperature_K is
        given, the rate is multiplied by a colour term: the ratio response / Bessell B photons of a
        blackbody at that temperature over the same ratio at T_VEGA_K (so temperature_K = T_VEGA_K
        changes nothing, and the Vega spectrum keeps its real shape, e.g. the Balmer jump).
        """
        wl = self.default_wavelengths.astype(float)
        resp = self.response(wl)
        phi = vega_photon_density(wl) * 10 ** (0.4 * VEGA_B_MAGNITUDE)
        rate = self.mirror_area_per_telescope[self.telescope_type] * 1e4 * \
            _trapz(resp * phi, wl) * 10    # 1e4: m2 -> cm2, 10: nm -> A
        if temperature_K is None:
            return rate
        s_b = filter_curve("bessell_B", wl) / wl  # photon-counting

        def response_over_b(T):
            bb = planck_photon(wl * nm_m, T)
            return _trapz(resp * bb, wl) / _trapz(s_b * bb, wl)

        return rate * response_over_b(temperature_K) / response_over_b(T_VEGA_K)

    def num_photons_from_mag(self, obs_times, mag, temperature_K=None):
        """
        Star photo-electrons in obs_times [ms] for a star of Johnson B magnitude `mag` (Vega system),
        whatever the optical filter. See star_photon_rate.
        """
        return obs_times * 1e-3 * self.star_photon_rate(temperature_K) * 10 ** (-0.4 * np.asarray(mag))

    @staticmethod
    def nsb_spectrum(moon=False):
        # Read nsb spectrum
        filename = "Spectra_NSB_ref.txt" if not moon else "Spectra_NSB_halfmoon.txt"
        data = ascii.read("{}/{}".format(reference_files_path, filename), format='no_header',
                        names = ['wavelength(nm)', 'flux'], data_start=0,
                        guess=False, delimiter=' ')
        # Note flux is in units of 1 / (nm m^2 ns sr)
        return data

    def pixel_solid_angle(self):
        """Solid angle of one pixel (deg^2): the FoV is its flat-to-flat size; hexagon, square or circle."""
        d = self.pixel_FoV_per_telescope[self.telescope_type]
        shape = self.pixel_shape_per_telescope.get(self.telescope_type, 'circle')
        area = {'hexagon': np.sqrt(3) / 2 * d ** 2, 'square': d ** 2}.get(shape, np.pi / 4 * d ** 2)
        return area * u.deg**2

    def psf_containment(self, psf_mm, convention):
        """
        Fraction of the light of a centred point source that falls inside one pixel, for a Gaussian PSF
        given by the logged number `psf_mm` (camera plane) under a stated convention:

          'r39' : 39 % containment RADIUS (a Gaussian holds 39.3 % within r = sigma): sigma = psf (MAGIC elog)
          'd80' : 80 % containment diameter (LST D80): sigma = psf / (2 sqrt(2 ln 5)) (LST-1 elog)

        An upper limit: mispointing lowers it. Returns 1.0 for psf_mm None/NaN.
        """
        if psf_mm is None or not np.isfinite(psf_mm):
            return 1.0
        sigmas = {'r39': psf_mm, 'd80': psf_mm / (2.0 * R80_OVER_SIGMA)}
        if convention not in sigmas:
            raise ValueError(f"unknown PSF convention {convention!r}; use 'r39' or 'd80'")
        sigma = sigmas[convention]
        flat_mm = self.pixel_FoV_per_telescope[self.telescope_type] * np.deg2rad(1.0) * \
            self.focal_length_mm_per_telescope[self.telescope_type]
        return pixel_containment(sigma, flat_mm, self.pixel_shape_per_telescope.get(self.telescope_type, 'hexagon'))

    def number_nsb_photons(self, moon=False):
        """
        Number of night sky background photons detected in 1ns

        Parameters:
        ----------
        moon: Boolean
        The default is set to no moon NSB, set to True for half moon observations

        """
        if not moon:
            nsb_table = self.nsb_spectrum(False)
        else:
            nsb_table = self.nsb_spectrum(True)
        if self.filter_type is None:
            transmission = self.total_telescope_transmission(nsb_table['wavelength(nm)'])
        else:
            transmission = self.optical_filter_transmission(nsb_table['wavelength(nm)']) * \
                            self.total_telescope_transmission(nsb_table['wavelength(nm)'])
        
        number_of_photons = np.sum(nsb_table['flux']*transmission)
        pixel_area = self.pixel_solid_angle()
        collection_area = self.mirror_area_per_telescope[self.telescope_type]*pixel_area.to('sr').value
        return number_of_photons*collection_area

    def signal_to_noise_ratio(self, magnitude, time_binning, ntels = 1, temperature_K=None):
        # Simple ratio of the photons detected from the star to the photons from the
        # night sky background (star of Johnson B magnitude `magnitude`) for the time binning used.
        # Uniform sampling required!
        # relative_intensities is the y axis of the diffraction pattern, relative to 'magnitude'
        #
        # This is the pure Poisson shot-noise (star + NSB) baseline. For bright
        # stars its fractional noise shrinks as ~1/sqrt(N_star) and the SNR is
        # wildly optimistic. See realistic_signal_to_noise_ratio / noise_budget
        # for the model that also includes the excess-noise factor, electronic
        # noise and atmospheric scintillation.
        if not isinstance(time_binning, u.Quantity):
            raise Exception('Wrong input', 'Times needs to have units')
        num_star_photons = self.num_photons_from_mag(time_binning.to('ms').value, magnitude, temperature_K)
        nsb_photons = self.number_nsb_photons(moon=False)*time_binning.to('ns').value
        # Noise level is shot noise (NSB+STAR), divided by sqrt(n telescopes)
        noise_level = np.sqrt(num_star_photons+nsb_photons)/np.sqrt(ntels)
        return num_star_photons/noise_level

    def noise_budget(self, magnitude, time_binning, ntels=1, airmass=1.0,
                     excess_noise_factor=None, pedestal_pe_per_ns=None,
                     aperture_diameter_m=None, scint_coeff=1.5, wind_speed_ms=20.0,
                     variance_floor=True, variance_k=None, n_frame_samples=None,
                     red_noise_rms=0.0, scintillation=True, moon=False, extinction=False,
                     throughput_factor=1.0, temperature_K=None, psf_containment=1.0):
        """
        Per-sample fractional-flux noise budget for a star of apparent
        `magnitude` observed at cadence `time_binning`.

        Unlike signal_to_noise_ratio (pure star+NSB shot noise), this combines
        physically-motivated terms in quadrature, all expressed as a
        fraction of the star signal S (so they are directly comparable to the
        measured std/mean of a real light curve):

          - shot noise, inflated by the photodetector excess-noise factor F:
                F**2 * (S + B) / S**2
          - atmospheric scintillation (Young 1967 / Osborn et al. 2015),
            magnitude-independent, see scintillation_variance: long exposure
                scint_coeff**2 * 1e-5 * D_m**(-4/3) * X**3 * exp(-2 h_obs / H) / t_s
            (D in metres), saturating at a t-independent value for t << D / v;

          - variance-estimator floor (the fast signal is a std^2 estimate from
            N samples per frame): 2 k / N per frame * ((S + B + P) / S)**2, with P
            the electronic pedestal (the estimator's relative error is relative to the
            whole variance, not to the star). Independent between
            frames, so it averages down with the number of frames in a bin
            (time_binning / frame duration) as well as with `ntels`;
          - red (correlated) noise of fractional rms `red_noise_rms` in the
            search band: not reduced by binning or by combining telescopes
            (the same wind/tracking/transparency hits all of them).

        Every white variance term is divided by `ntels`, matching
        signal_to_noise_ratio. Note scintillation is only partially independent
        between the ~100 m-separated MAGIC/LST dishes, so 1/ntels is a best
        case for that term.

        Parameters
        ----------
        magnitude : float
            Johnson B magnitude of the star (Vega system), whatever `filter_type` is.
        time_binning : astropy.units.Quantity
            Sampling cadence (needs time units).
        ntels : int
            Number of telescopes combined.
        airmass : float
            Airmass X = sec(z) of the target. Default 1.0 (zenith).
        excess_noise_factor, pedestal_pe_per_ns, aperture_diameter_m : float, optional
            Override the per-telescope defaults (excess_noise_factor_per_telescope,
            pedestal_pe_per_ns_per_telescope, and 2*sqrt(area/pi) respectively).
        scint_coeff : float
            Empirical scaling on the scintillation sigma (Young's coefficient
            is for a median site; Osborn et al. 2015 find ~1.5 at La Palma,
            the default).
        wind_speed_ms : float
            Turbulence-weighted wind speed (assumed) that sets the short/long-exposure
            crossover t ~ D / v of the scintillation.
        variance_floor : bool
            Include the variance-estimator floor (default True; zero for
            telescopes with k = 0).
        variance_k, n_frame_samples : float, optional
            Override k and N of the floor.
        red_noise_rms : float
            Fractional rms of the correlated noise in the signal band
            (e.g. a few 1e-3 after a 3 Hz high-pass; see inband_rms).
        scintillation : bool
            Set False to drop the scintillation term entirely.
        moon : bool
            Passed to number_nsb_photons for the NSB level.
        extinction : bool
            If True, also evaluate the atmospheric extinction of the star light
            at `airmass` (exp(-tau * X), see set_airmass). By default `airmass`
            only enters the scintillation term and the star photons are
            computed at the zenith, as before.
        throughput_factor : float
            Extra multiplicative efficiency applied to both the star and the
            NSB photons, e.g. the fraction of mirror area left when some
            panels are faulty, or a measured cloud transmission.
        temperature_K : float, optional
            Effective temperature of the star, for the colour term of star_photon_rate
            (default: Vega colours, B - V = 0). See physics.teff_from_spectral_type.
        psf_containment : float
            Fraction of the star light inside the pixel (see Instrument.psf_containment).
            It multiplies the star photons only: the NSB is diffuse, so the light the PSF
            spreads out of the pixel is replaced by the sky spread in from its neighbours.
            (throughput_factor, which is collecting area, scales both.)

        Returns
        -------
        dict with keys star_photons (photo-electrons: the response includes the QE), nsb_photons, sigma_shot,
        sigma_scint and sigma_total (all sigmas in fractional-flux units).
        """
        if not isinstance(time_binning, u.Quantity):
            raise Exception('Wrong input', 'Times needs to have units')

        if excess_noise_factor is None:
            excess_noise_factor = self.excess_noise_factor_per_telescope[self.telescope_type]
        if pedestal_pe_per_ns is None:
            pedestal_pe_per_ns = self.pedestal_pe_per_ns_per_telescope[self.telescope_type]
        if aperture_diameter_m is None:
            area = self.mirror_area_per_telescope[self.telescope_type]
            aperture_diameter_m = 2.0 * np.sqrt(area / np.pi)

        previous_airmass = self.airmass
        if extinction:
            self.set_airmass(airmass)
        try:
            S = self.num_photons_from_mag(time_binning.to('ms').value, magnitude, temperature_K)
        finally:
            self.set_airmass(previous_airmass)
        B = self.number_nsb_photons(moon=moon) * time_binning.to('ns').value
        P = pedestal_pe_per_ns * time_binning.to('ns').value * throughput_factor
        S = S * throughput_factor * psf_containment
        B = B * throughput_factor

        # Fractional variances (relative to the star signal S), combined
        # across telescopes as 1/ntels.
        var_shot = excess_noise_factor**2 * (S + B) / S**2 / ntels

        if scintillation:
            var_scint = scintillation_variance(
                aperture_diameter_m, time_binning.to('s').value, airmass, self.site_altitude_m,
                scint_coeff, wind_speed_ms) / ntels
        else:
            var_scint = 0.0

        var_floor = 0.0
        if variance_floor:
            k = self.variance_estimator_k_per_telescope[self.telescope_type] if variance_k is None else variance_k
            n = self.variance_estimator_n_samples if n_frame_samples is None else n_frame_samples
            frame_s = FRAME_SAMPLES_DURATION_S
            n_frames = max(1.0, time_binning.to('s').value / frame_s)
            var_floor = 2.0 * k / n / n_frames / ntels * ((S + B + P) / S) ** 2
        var_red = red_noise_rms ** 2

        sigma_total = np.sqrt(var_shot + var_scint + var_floor + var_red)

        return {
            'star_photons': S,
            'nsb_photons': B,
            'sigma_shot': np.sqrt(var_shot),
            'sigma_scint': np.sqrt(var_scint),
            'sigma_floor': np.sqrt(var_floor),
            'sigma_red': np.sqrt(var_red),
            'sigma_total': sigma_total,
        }

    def realistic_signal_to_noise_ratio(self, magnitude, time_binning, ntels=1, **kwargs):
        """
        Signal-to-noise ratio using the full noise_budget (star+NSB shot noise
        with excess-noise factor, electronic noise and atmospheric
        scintillation), i.e. SNR = 1 / sigma_total.

        Same required arguments as signal_to_noise_ratio; extra keyword
        arguments (airmass, scint_coeff, excess_noise_factor, pedestal_pe_per_ns, psf_containment,
        aperture_diameter_m, scintillation, moon, extinction, throughput_factor,
        temperature_K)
        are forwarded to noise_budget.
        """
        budget = self.noise_budget(magnitude, time_binning, ntels=ntels, **kwargs)
        return 1.0 / budget['sigma_total']