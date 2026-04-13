import numpy as np


def load_response_file(filename):
    """
    Load wavelength-dependent response from a text file.

    Expected format:
    wavelength_nm   transmission_percent
    """

    data = np.loadtxt(filename)

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