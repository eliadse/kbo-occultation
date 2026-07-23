# examples/preprocess_observations.py
"""
Pre-process the raw stat-binary observations into compact ``.npz`` caches.

Each raw ``.bin`` is parsed once; only the information the package actually
uses (the per-channel variance flux proxy and the start time) is kept, the
outlier correction is applied, and both the raw and cleaned variance are stored
in a small ``.npz`` beside the source file. Downstream workflows then load the
cache (``LightCurve.from_preprocessed*`` / ``search_observation`` on the
``.npz``) instead of re-parsing the ~140 MB binary every run.

Re-run this after new data lands, or after changing the outlier parameters.
"""

import glob
import os

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.photometry import preprocess_observation

DATA_DIR = f"{PACKAGE_DATA}/observations"


def main():
    bin_files = sorted(glob.glob(f"{DATA_DIR}/Spectrum*.bin"))
    if not bin_files:
        print(f"No .bin files found under {DATA_DIR}")
        return

    for bin_path in bin_files:
        out_path = preprocess_observation(bin_path)
        raw_mb = os.path.getsize(bin_path) / 1e6
        npz_mb = os.path.getsize(out_path) / 1e6
        ratio = raw_mb / npz_mb if npz_mb else float("nan")
        print(f"{os.path.basename(bin_path)}: "
              f"{raw_mb:.1f} MB -> {npz_mb:.1f} MB ({ratio:.1f}x) "
              f"-> {os.path.basename(out_path)}")


if __name__ == "__main__":
    main()
