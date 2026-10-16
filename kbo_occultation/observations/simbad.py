"""
Star coordinates, magnitudes and spectral types from SIMBAD, for filling the
observation log automatically (``kbo-night prepare``). Needs network access
and astroquery; both are only used here, at prepare time.
"""

import re
import warnings
from typing import Dict, Iterable, Optional

# HD-4335, HD170263, HIP_12345 -> catalogue + number
_CATALOGUE_NUMBER = re.compile(r"^([A-Za-z]+)[-_ ]*(\d.*)$")


def simbad_name(name: str) -> str:
    """SIMBAD identifier for a run name from a filename: ``HD-4335`` / ``HD170263`` -> ``HD 4335`` / ``HD 170263``."""
    m = _CATALOGUE_NUMBER.match(name.strip())
    if m:
        return f"{m.group(1)} {m.group(2)}"
    return re.sub(r"[-_]+", " ", name).strip()


def _value(row, col, cast):
    import numpy as np

    v = row[col]
    if np.ma.is_masked(v):
        return None
    v = cast(v)
    if isinstance(v, str):
        v = v.strip()
        return v or None
    return None if v != v else v   # NaN -> None


def query_stars(names: Iterable[str]) -> Dict[str, Optional[dict]]:
    """
    ``{name: {ra_deg, dec_deg, mag_B, mag_V, spectral_type, simbad_id}}`` for each
    run name (resolved with ``simbad_name``); None for names SIMBAD does not
    know. If SIMBAD cannot be reached, warns and returns None for every name.
    """
    from astroquery.simbad import Simbad

    names = list(names)
    if not names:
        return {}
    simbad = Simbad()
    simbad.add_votable_fields("B", "V", "sp_type")
    try:
        table = simbad.query_objects([simbad_name(n) for n in names])
    except Exception as e:  # network / service errors: leave the fields for the next run
        warnings.warn(f"SIMBAD query failed ({e}); coordinates and magnitudes left unset")
        return {n: None for n in names}

    out = {n: None for n in names}
    for i, row in enumerate(table):
        # object_number_id is the 1-based position of the queried name
        idx = int(row["object_number_id"]) - 1 if "object_number_id" in table.colnames else i
        name = names[idx]
        ra = _value(row, "ra", float)
        if ra is None:
            continue
        out[name] = {
            "ra_deg": round(ra, 5),
            "dec_deg": round(_value(row, "dec", float), 5),
            "mag_B": _round(_value(row, "B", float)),
            "mag_V": _round(_value(row, "V", float)),
            "spectral_type": _value(row, "sp_type", str),
            "simbad_id": " ".join(str(row["main_id"]).split()),
        }
    return out


def _round(v, ndigits=3):
    return None if v is None else round(v, ndigits)
