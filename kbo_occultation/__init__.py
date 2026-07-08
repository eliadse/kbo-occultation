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
)

from .sweep import (
    SweepPoint,
    run_parameter_sweep,
    to_records,
    to_dataframe,
)

from .detectability import (
    spatial_to_time,
    resample_to_cadence,
    peak_snr,
    matched_filter_snr,
    event_duration,
    is_detectable,
)
