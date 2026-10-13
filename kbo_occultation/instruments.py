"""
Functions and reference files created and compiled by Tarek Hassan
"""
import numpy as np
from astropy.io import ascii
from astropy import units as u
from scipy import interpolate

from kbo_occultation import PACKAGE_DATA
reference_files_path = "{}/iact_reference_values".format(PACKAGE_DATA)

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

#TODO Add a function that lists all available filters and/or transmission files
class Instrument():
    """

    A class to compute the transmission of the telescopes and their signal to noise ratio

    """
    def __init__(self, telescope_type='MAGIC', filter_type=None, site='Paranal'):
        self.telescope_type = telescope_type
        self.filter_type = filter_type
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

        nsb_per_telescope: Night-sky background: [photons/ns/pixel]
        
        pixel_FoV_per_telescope: Telescope pixel FoV (diameter)[deg]
        
        mirror_reflectivity_per_telescope: Telescope shadowing?
        
        qes_per_telescope: Quantum efficiency of the telescopes
        
        transmission_per_telescope: Camera elements transmission (plexiglass, etc.)
        
        funnel_transmission_per_telescope: Funnel efficiency at 0 offset
        
        camera_transmission_per_telescope: Camera transmission of the telescopes
        
        mirror_area_per_telescope: Geometric mirror area BEFORE SHADOWING [m^2]
        
        number_of_photons_from_Vega: In different bands [photons cm-2 s-1 A-1]

        """
        self.default_wavelengths = np.arange(100, 1000, 1)

        self.nsb_per_telescope = {
            'MAGIC': 0.09770744068777938,
            'LST': 0.26516,
            'MST-N': 0.25612,
            'pSCT': 0.04711,
            'SST-1M': 0.04068,
            'VERITAS': 0.0662455616546292,
            'VERITAS++': 0.0662455616546292
            }

        # LST NSB photon rate is 0.26516458670986737 ph/ns/pixel
        # MST-N NSB photon rate is 0.25612364288761486 ph/ns/pixel
        # pSCT NSB photon rate is 0.017265456577839695 ph/ns/pixel
        # SST-1M NSB photon rate is 0.04004496264110617 ph/ns/pixel
        # VERITAS NSB photon rate is 0.08898415585269348 ph/ns/pixel

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
            'LST': 0.892,        #tweaked to match the NSB rate used
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
        
        # We extract the number of photons from Vega: http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
        # [photons cm-2 s-1 A-1]
        self.number_of_photons_from_Vega = {
            'u': 1539.3,
            'g': 1134.6,
            'r': 875.4,
            'i': 714.5,
            'z': 602.2,
            'MAGIC_SII': 356.0,
            'QE': 2020.2,
            'MAGIC_QE': 2020.2,
            'gaia_bp': 9.73
            }

        # Photodetector excess-noise factor F (variance is inflated by F^2
        # relative to pure Poisson shot noise). Values from MAGIC-LST.yml.
        self.excess_noise_factor_per_telescope = {
            'MAGIC': 1.15,
            'LST': 1.15,
            'MST-N': 1.15,
            'pSCT': 1.15,
            'SST-1M': 1.15,
            'VERITAS': 1.15,
            'VERITAS++': 1.15
            }

        # Electronic/readout noise, RMS photoelectrons per sample. From MAGIC-LST.yml.
        self.elec_noise_per_telescope = {
            'MAGIC': 0.7,
            'LST': 0.7,
            'MST-N': 0.7,
            'pSCT': 0.7,
            'SST-1M': 0.7,
            'VERITAS': 0.7,
            'VERITAS++': 0.7
            }

        # Observatory altitude (m) for the atmospheric scintillation term.
        # Paranal ~2635 m, ORM/La Palma ~2200 m.
        self.site_altitude_m = 2635.0 if self.site == 'Paranal' else 2200.0

        self.transmission = self.total_telescope_transmission()

    def print_parameters(self):
        print(f"Telescope type: {self.telescope_type}")
        print(f"Optical filter: {self.filter_type}")
        print(f"Telescope site: {self.site}")

    def set_filter(self, filter_type):
        """
        Set or change the optical filter of the telescope.
        """
        self.filter_type = filter_type

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
        result[mask] = np.exp(-1 * f(wavelength[mask]))
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
        # Interpolate optical filter transmission for any filter available
        # The guess=True addition let's us use any filter file
        data = ascii.read("{}/optical_filters/optical_filter_{}.txt".format(PACKAGE_DATA, self.filter_type),
                        guess=True)
        f = interpolate.interp1d(data['col1'], data['col2'])
        result = np.zeros_like(wavelength, dtype=float)
        mask = np.logical_and(wavelength >= np.min(data['col1']), wavelength <= np.max(data['col1']))
        result[mask] = f(wavelength[mask])
        return result/100.

    def response_band(self, threshold=0.01, margin_nm=10.0, n_lambda=60):
        """
        Auto-detect the useful wavelength range for an instrument+filter
        combo (where its total throughput is above `threshold` of its peak),
        instead of guessing a bandpass by hand.
        """
        wl = self.default_wavelengths.astype(float)
        resp = (self.total_transmission_filter(wl) if self.filter_type
                else self.total_transmission(wl))
        peak = resp.max()
        above = wl[resp > threshold * peak]
        lam_min = max(wl.min(), above.min() - margin_nm)
        lam_max = min(wl.max(), above.max() + margin_nm)
        return float(lam_min), float(lam_max), n_lambda

# TODO Move the Vega stuff to physics instead. It doesn't make sense to calculate the photons
#from Vega just using a filter and not the full instrument transmission, so fix that too 
    @staticmethod
    def star_flux(wavelength, T):
        # These are not needed anymore, because I directly calculate everything from Vega filters
        h = 6.626e-34
        c = 3.0e+8
        k = 1.38e-23
        a = 2.0*h*c**2
        wavelength = wavelength*1e-9
        b = h*c/(wavelength*k*T)
        intensity = a/ ( (wavelength**5) * (np.exp(b) - 1.0) )
        return intensity

    def star_flux_folded_with_transmission(self, wavelength, T):
        # These are not needed anymore, because I directly calculate everything from Vega filters
        selected_filter = 1. if self.filter_type is None else self.optical_filter_transmission(wavelength)
        return star_flux(wavelength, T) * self.total_transmission(wavelength) * selected_filter

    def num_photons_from_mag(self, obs_times, mag):
        """
            The number of photons for a given magnitude is interpolated from the number of photons of Vega.
            Functions are calculated from 100 to 1000nm
        """
        if self.filter_type is None:
            filter_values = self.total_transmission()
            photons_from_Vega = self.number_of_photons_from_Vega['QE']

        else:
            filter_values = self.optical_filter_transmission() * self.total_transmission()
            photons_from_Vega = self.number_of_photons_from_Vega[self.filter_type]

        f = interpolate.interp1d(self.default_wavelengths, filter_values)
        #integral = integrate.quad(f, np.min(self.default_wavelengths), np.max(self.default_wavelengths), limit=60)
        # The quad integration method struggles with functions that are flat and have steep
        # edges, so it's better to switch to trapz
        trapz = getattr(np, "trapezoid", None) or np.trapz  # np.trapz removed in NumPy 2.0
        integral = trapz(f.y, f.x)
        # 1e3 because going from ms to s
        # 1e4 because going from m2 to cm2
        # 10  because going from nm to A
        num_photons_over_filter = obs_times * 1e-3 * \
                                self.mirror_area_per_telescope[self.telescope_type] * 1e4 * \
                                integral * 10 * \
                                photons_from_Vega

        return num_photons_over_filter / (2.512 ** mag)

    @staticmethod
    def nsb_spectrum(moon=False):
        # Read nsb spectrum
        filename = "Spectra_NSB_ref.txt" if not moon else "Spectra_NSB_halfmoon.txt"
        data = ascii.read("{}/{}".format(reference_files_path, filename), format='no_header',
                        names = ['wavelength(nm)', 'flux'], data_start=0,
                        guess=False, delimiter=' ')
        # Note flux is in units of 1 / (nm m^2 ns sr)
        return data

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
        pixel_area = np.pi*np.power(self.pixel_FoV_per_telescope[self.telescope_type]/2., 2) * u.deg**2
        collection_area = self.mirror_area_per_telescope[self.telescope_type]*pixel_area.to('sr').value
        return number_of_photons*collection_area

    def signal_to_noise_ratio(self, magnitude, time_binning, ntels = 1):
        # Simple ratio of the photons detected from the star to the photons from the
        # night sky background magnitude at a given filter for the time binning used.
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
        num_star_photons = self.num_photons_from_mag(time_binning.to('ms').value, magnitude)
        nsb_photons = self.number_nsb_photons(moon=False)*time_binning.to('ns').value
        # Noise level is shot noise (NSB+STAR), divided by sqrt(n telescopes)
        noise_level = np.sqrt(num_star_photons+nsb_photons)/np.sqrt(ntels)
        return num_star_photons/noise_level

    def noise_budget(self, magnitude, time_binning, ntels=1, airmass=1.0,
                     excess_noise_factor=None, elec_noise_pe=None,
                     aperture_diameter_m=None, scint_coeff=1.0,
                     scintillation=True, moon=False):
        """
        Per-sample fractional-flux noise budget for a star of apparent
        `magnitude` observed at cadence `time_binning`.

        Unlike signal_to_noise_ratio (pure star+NSB shot noise), this combines
        three physically-motivated terms in quadrature, all expressed as a
        fraction of the star signal S (so they are directly comparable to the
        measured std/mean of a real light curve):

          - shot noise, inflated by the photodetector excess-noise factor F:
                F**2 * (S + B) / S**2
          - electronic/readout noise of `elec_noise_pe` RMS photoelectrons:
                elec_noise_pe**2 / S**2
          - atmospheric scintillation (Young 1967), magnitude-independent:
                scint_coeff**2 * 10e-6 * D_cm**(-4/3) * X**3
                                * exp(-2 h_obs / H) / (2 t_s)
            with D the aperture diameter, X the airmass, h_obs the site
            altitude and H = 8000 m the atmospheric scale height.

        Every variance term is divided by `ntels`, matching
        signal_to_noise_ratio. Note scintillation is only partially independent
        between the ~100 m-separated MAGIC/LST dishes, so 1/ntels is a best
        case for that term.

        Parameters
        ----------
        magnitude : float
            Apparent magnitude in the current filter (unfiltered if None).
        time_binning : astropy.units.Quantity
            Sampling cadence (needs time units).
        ntels : int
            Number of telescopes combined.
        airmass : float
            Airmass X = sec(z) of the target. Default 1.0 (zenith).
        excess_noise_factor, elec_noise_pe, aperture_diameter_m : float, optional
            Override the per-telescope defaults (excess_noise_factor_per_telescope,
            elec_noise_per_telescope, and 2*sqrt(area/pi) respectively).
        scint_coeff : float
            Empirical scaling on the scintillation sigma (Young's coefficient
            is for a median site; Osborn et al. 2015 find ~1.5 at La Palma).
        scintillation : bool
            Set False to drop the scintillation term entirely.
        moon : bool
            Passed to number_nsb_photons for the NSB level.

        Returns
        -------
        dict with keys star_photons, nsb_photons, sigma_shot, sigma_elec,
        sigma_scint and sigma_total (all sigmas in fractional-flux units).
        """
        if not isinstance(time_binning, u.Quantity):
            raise Exception('Wrong input', 'Times needs to have units')

        if excess_noise_factor is None:
            excess_noise_factor = self.excess_noise_factor_per_telescope[self.telescope_type]
        if elec_noise_pe is None:
            elec_noise_pe = self.elec_noise_per_telescope[self.telescope_type]
        if aperture_diameter_m is None:
            area = self.mirror_area_per_telescope[self.telescope_type]
            aperture_diameter_m = 2.0 * np.sqrt(area / np.pi)

        S = self.num_photons_from_mag(time_binning.to('ms').value, magnitude)
        B = self.number_nsb_photons(moon=moon) * time_binning.to('ns').value

        # Fractional variances (relative to the star signal S), combined
        # across telescopes as 1/ntels.
        var_shot = excess_noise_factor**2 * (S + B) / S**2 / ntels
        var_elec = elec_noise_pe**2 / S**2 / ntels

        if scintillation:
            t_s = time_binning.to('s').value
            D_cm = aperture_diameter_m * 100.0
            H = 8000.0  # atmospheric scale height (m)
            var_scint = (scint_coeff**2 * 10e-6 * D_cm**(-4.0 / 3.0) *
                         airmass**3 * np.exp(-2.0 * self.site_altitude_m / H) /
                         (2.0 * t_s)) / ntels
        else:
            var_scint = 0.0

        sigma_total = np.sqrt(var_shot + var_elec + var_scint)

        return {
            'star_photons': S,
            'nsb_photons': B,
            'sigma_shot': np.sqrt(var_shot),
            'sigma_elec': np.sqrt(var_elec),
            'sigma_scint': np.sqrt(var_scint),
            'sigma_total': sigma_total,
        }

    def realistic_signal_to_noise_ratio(self, magnitude, time_binning, ntels=1, **kwargs):
        """
        Signal-to-noise ratio using the full noise_budget (star+NSB shot noise
        with excess-noise factor, electronic noise and atmospheric
        scintillation), i.e. SNR = 1 / sigma_total.

        Same required arguments as signal_to_noise_ratio; extra keyword
        arguments (airmass, scint_coeff, excess_noise_factor, elec_noise_pe,
        aperture_diameter_m, scintillation, moon) are forwarded to noise_budget.
        """
        budget = self.noise_budget(magnitude, time_binning, ntels=ntels, **kwargs)
        return 1.0 / budget['sigma_total']