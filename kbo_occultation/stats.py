# kbo_occultation/stats.py
"""
Small statistics helpers shared by the observations and analysis
subpackages.
"""

import numpy as np


def robust_sigma(x: np.ndarray) -> float:
    """
    MAD-based noise sigma, robust to a brief injected/real dip (or a
    handful of outliers) not biasing the estimate the way a plain
    ``std`` would over a long series with a short event.
    """
    x = np.asarray(x)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))
