# kbo_occultation/config.py

from dataclasses import dataclass

@dataclass
class KBOConfig:
    radius_m: float
    distance_au: float
    impact_parameter_m: float = 0.0

@dataclass
class StarConfig:
    temperature_K: float
    angular_radius_mas: float

@dataclass
class BandpassConfig:
    lam_min_nm: float
    lam_max_nm: float
    n_lambda: int

@dataclass
class GridConfig:
    x_max_m: float
    n_x: int

@dataclass
class NumericalConfig:
    n_int: int = 40
    n_r_grid: int = 3000
    n_star_side: int = 32

standard_sampling = 2**18
standard_sample_duration = 2 #ns
lightcurve_sampling = 50 # [Hz]. 2 * Nyquist freq. From plot_diffraction_psd with r_* = 0.03mas
