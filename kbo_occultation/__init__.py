import os

PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))
PACKAGE_DATA = os.path.join(PACKAGE_ROOT, 'data')

from .simulation import (
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

from .sweep import (
    SweepPoint,
    run_parameter_sweep,
    to_records,
    to_dataframe,
)

from .detectability import (
    spatial_to_time,
    spatial_power_spectrum,
    nyquist_frequency,
    resample_to_cadence,
    bin_average,
    peak_snr,
    matched_filter_snr,
    event_duration,
    is_detectable,
)

from .injection import (
    InjectionResult,
    inject_occultation,
    random_injection_time,
)

from .matched_filter import (
    MatchedFilterResult,
    Candidate,
    robust_sigma,
    compute_baseline_and_sigma,
    sliding_matched_filter_snr,
    shape_veto_chi2,
    find_candidates,
)

from .injection_batch import (
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

from .kinematics import (
    ShadowKinematics,
    shadow_velocity_opposition,
    shadow_velocity_from_target,
)

from .coincidence import (
    TelescopeSearchResult,
    CoincidenceEvent,
    AccidentalRateResult,
    coincidence_tolerance_s,
    match_coincidences,
    accidental_coincidence_rate,
)

from .upper_limits import (
    EfficiencyCurve,
    SurveyExposure,
    UpperLimitResult,
    fresnel_scale_m,
    effective_shadow_width_m,
    poisson_upper_limit,
    exposure_omega_deg2,
    cumulative_density_upper_limit,
)

from .search import (
    ObservationSearchResult,
    search_observation,
)

from .dc_combine import (
    load_dc_series,
    interpolate_dc_to_fast,
    combine_fast_and_dc,
    combine_lightcurve_with_dc,
    smooth_dc_series,
    detrend_lightcurve_with_dc,
    flag_dc_excursions,
    despike_lightcurve_with_dc,
)

from .filtering import (
    highpass_fft,
    highpass_butterworth,
    highpass_lightcurve,
    remove_outliers_lightcurve,
)
