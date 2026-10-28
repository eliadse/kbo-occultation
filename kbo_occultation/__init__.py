import os

PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))
PACKAGE_DATA = os.path.join(PACKAGE_ROOT, 'data')

from .simulation.engine import (
    compute_lightcurve,
    simulate_poly_point,
    OccultationEngine,
)
#from .instruments import (
#    load_response_file,
#    build_response_function,
#    combine_responses
#)

from .config import (
    KBOConfig,
    StarConfig,
    BandpassConfig,
    GridConfig,
    NumericalConfig,
    TelescopeConfig,
    ArrayConfig,
    ORM_ARRAY,
)

from .simulation.sweep import (
    SweepPoint,
    run_parameter_sweep,
    to_records,
    to_dataframe,
)

from .simulation.detectability import (
    spatial_to_time,
    spatial_power_spectrum,
    cumulative_quantile_frequency,
    nyquist_frequency,
    resample_to_cadence,
    bin_average,
    peak_snr,
    matched_filter_snr,
    event_duration,
    is_detectable,
)

from .analysis.injection import (
    InjectionResult,
    inject_occultation,
    random_injection_time,
)

from .analysis.matched_filter import (
    MatchedFilterResult,
    Candidate,
    robust_sigma,
    compute_baseline_and_sigma,
    sliding_matched_filter_snr,
    shape_veto_chi2,
    calibrate_chi2_threshold,
    empirical_filter_sigma,
    find_candidates,
)

from .analysis.injection_batch import (
    TrialResult,
    MethodNoiseProfile,
    ArrayTrialResult,
    build_search_template,
    run_injection_monte_carlo,
    run_array_injection_monte_carlo,
    efficiency_curve,
    to_records as batch_to_records,
    to_dataframe as batch_to_dataframe,
)

from .simulation.kinematics import (
    ShadowKinematics,
    shadow_velocity_opposition,
    shadow_velocity_from_target,
)

from .analysis.coincidence import (
    TelescopeSearchResult,
    CoincidenceEvent,
    AccidentalRateResult,
    coincidence_tolerance_s,
    match_coincidences,
    accidental_coincidence_rate,
)

from .analysis.upper_limits import (
    EfficiencyCurve,
    SurveyExposure,
    UpperLimitResult,
    fresnel_scale_m,
    effective_shadow_width_m,
    poisson_upper_limit,
    exposure_omega_deg2,
    cumulative_density_upper_limit,
)

from .analysis.rates import (
    event_rate,
    kbo_velocity,
    surface_density,
    mean_diameter,
    effective_shadow_width,
    Q_1,
    Q_L,
)

from .analysis.search import (
    ObservationSearchResult,
    search_observation,
)

from .observations.dc_combine import (
    load_dc_series,
    interpolate_dc_to_fast,
    combine_fast_and_dc,
    combine_lightcurve_with_dc,
    smooth_dc_series,
    detrend_lightcurve_with_dc,
    flag_dc_excursions,
    despike_lightcurve_with_dc,
)

from .observations.filtering import (
    highpass_fft,
    highpass_butterworth,
    highpass_lightcurve,
    remove_outliers_lightcurve,
)

from .observations.photometry import (
    LightCurve,
    preprocess_observation,
    good_time_lightcurve,
)

from .observations.log import (
    Observation,
    Night,
    ObservationLog,
    load_observation_log,
    validate,
)

from .analysis.snr import snr_report
