"""
Observing conditions for the real observations: parsing of the MAGIC central
control (CC) ``.rep`` reports, and pointing/moon geometry computed from the
target coordinates.

The CC reports are ~300 MB per telescope and night and live outside the repo,
so ``kbo-night prepare`` (observations.prepare.extract_conditions) parses them
once with ``read_cc_reports`` and caches the observing window as small CSVs
(``save_conditions`` / ``load_conditions``). ``conditions_table`` turns those
caches into one row of conditions per observation.

Report fields used (positions checked against real 2025/2026 lines):

- ``CC-REPORT``: the ten floats before ``SCHEDULE`` are
  Zd, Az, target Dec, target RA, temperature, pressure, wind speed, humidity,
  (two unused); ``SCHEDULE`` is followed by the scheduled source. Zd/Az agree
  with the DRIVE-REPORT nominal position.
- ``WEATHER-REPORT``: ``CODE value`` pairs. TE temperature (C), DR pressure
  (hPa), FE humidity (%), WR wind direction (deg), WD average wind speed,
  WS peak wind speed, WG instantaneous wind speed (km/h). The CC-REPORT wind
  tracks WD. TK/TD/TS/RA/RS are kept with their raw codes.
- ``PYRO-REPORT`` (pyrometer): sky temperature (K), cloudiness (%), lid flag,
  air temperature (C) -- the air temperature matches WEATHER TE.
- ``AMC-REPORT``: the four integers after ``PANELS`` are kept raw as
  ``amc_panels_0..3`` = (enabled, undefined, adjusted OK, error) panels.
  ``amc_panels_2`` is only non-zero while an adjustment is running, and then
  ``amc_panels_2 + amc_panels_3 == amc_panels_0``; the error count stays
  constant between adjustments. Panels not correctly focused = undefined +
  error (``amc_unfocused_panels``).
- ``STARG-REPORT`` (starguider): Zd/Az deviation, number of identified stars
  and sky brightness (Mars MReportStarguider layout).
- ``LIDAR-REPORT``: only the ERROR_CODE is kept (the LIDAR was not
  operating during any of the campaign nights).
"""

import glob
import os
from typing import Dict, Optional

import numpy as np
import pandas as pd
from astropy.time import Time

from kbo_occultation import config

REPORT_PREFIXES = ("CC-REPORT", "WEATHER-REPORT", "PYRO-REPORT", "AMC-REPORT",
                   "STARG-REPORT", "LIDAR-REPORT")

WEATHER_CODES = {
    "TE": "temperature_C",
    "DR": "pressure_hPa",
    "FE": "humidity_pct",
    "WR": "wind_dir_deg",
    "WD": "wind_avg_kmh",
    "WS": "wind_peak_kmh",
    "WG": "wind_kmh",
}


def _report_time(tokens):
    """Unix time from the 7 date tokens (Y M D h m s ms) following the report status."""
    y, mo, d, h, mi, s, ms = tokens[3:10]
    return f"{y}-{mo}-{d}T{h}:{mi}:{s}.{ms}"


def _parse_cc(tokens):
    i = tokens.index("SCHEDULE")
    zd, az, dec, ra, temp, pressure, wind, hum = (float(x) for x in tokens[i - 10:i - 2])
    moon = tokens[tokens.index("GRB") - 1] if "GRB" in tokens else None
    if zd == 0.0 and az == 0.0:
        # CC writes 00.00 when it has no pointing to report (e.g. throughout
        # 2025-12-16); treat it as missing rather than as the zenith.
        zd = az = np.nan
    return {"zd_deg": zd, "az_deg": az, "target_dec_deg": dec, "target_ra_deg": ra,
            "temperature_C": temp, "pressure_hPa": pressure, "wind_kmh": wind,
            "humidity_pct": hum, "source": tokens[i + 1], "moon_flag": moon}


def _parse_weather(tokens):
    row = {}
    body = tokens[10:tokens.index("OVER")] if "OVER" in tokens else tokens[10:]
    for code, value in zip(body[::2], body[1::2]):
        row[WEATHER_CODES.get(code, f"weather_{code}")] = float(value)
    return row


def _parse_pyro(tokens):
    sky_k, cloud, lid, air = tokens[tokens.index("OVER") - 4:tokens.index("OVER")]
    return {"sky_temperature_K": float(sky_k), "cloudiness_pct": float(cloud),
            "lid_flag": int(lid), "air_temperature_C": float(air)}


def _parse_amc(tokens):
    i = tokens.index("PANELS")
    return {f"amc_panels_{k}": int(tokens[i + 1 + k]) for k in range(4)}


def _parse_starg(tokens):
    return {"starg_dzd_deg": float(tokens[18]), "starg_daz_deg": float(tokens[19]),
            "starg_n_stars": int(tokens[30]), "starg_sky_brightness": float(tokens[31])}


def _parse_lidar(tokens):
    code = tokens[tokens.index("ERROR_CODE") + 1] if "ERROR_CODE" in tokens else None
    return {"lidar_error_code": int(code) if code is not None else None}


_PARSERS = {
    "CC-REPORT": ("cc", _parse_cc),
    "WEATHER-REPORT": ("weather", _parse_weather),
    "PYRO-REPORT": ("pyro", _parse_pyro),
    "AMC-REPORT": ("amc", _parse_amc),
    "STARG-REPORT": ("starg", _parse_starg),
    "LIDAR-REPORT": ("lidar", _parse_lidar),
}


def read_cc_reports(rep_path: str, t_start: Optional[float] = None,
                    t_stop: Optional[float] = None) -> Dict[str, pd.DataFrame]:
    """
    Parse one CC ``.rep`` file into one DataFrame per report type
    (keys: cc, weather, pyro, amc, starg, lidar), each with a ``time`` column
    in unix seconds and the report's ``telescope`` (M0 = site-wide, M1, M2).

    ``t_start``/``t_stop`` (unix s) restrict the output to a time window; the
    file is still streamed once, but the huge DRIVE/CAMERA lines are skipped
    on their prefix without being split. Malformed lines are skipped.
    """
    rows = {name: [] for name, _ in _PARSERS.values()}
    with open(rep_path, errors="replace") as f:
        for line in f:
            if not line.startswith(REPORT_PREFIXES):
                continue
            tokens = line.split()
            name, parser = _PARSERS[tokens[0]]
            try:
                iso = _report_time(tokens)
                row = parser(tokens)
            except (ValueError, IndexError):
                continue
            row["iso"] = iso
            row["telescope"] = tokens[1]
            rows[name].append(row)

    out = {}
    for name, r in rows.items():
        df = pd.DataFrame(r)
        if df.empty:
            out[name] = df
            continue
        df.insert(0, "time", Time(list(df.pop("iso")), format="isot", scale="utc").unix)
        if t_start is not None:
            df = df[df["time"] >= t_start]
        if t_stop is not None:
            df = df[df["time"] <= t_stop]
        out[name] = df.reset_index(drop=True)
    return out


def save_conditions(reports: Dict[str, pd.DataFrame], out_dir: str, telescope: str):
    """Write ``read_cc_reports`` output as ``<out_dir>/<telescope>_<report>.csv``."""
    os.makedirs(out_dir, exist_ok=True)
    for name, df in reports.items():
        df.to_csv(os.path.join(out_dir, f"{telescope}_{name}.csv"), index=False)


def load_conditions(conditions_dir: str) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Inverse of ``save_conditions`` for every telescope: {telescope: {report: DataFrame}}."""
    out = {}
    for path in sorted(glob.glob(os.path.join(conditions_dir, "*_*.csv"))):
        telescope, name = os.path.splitext(os.path.basename(path))[0].split("_", 1)
        try:
            df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()
        out.setdefault(telescope, {})[name] = df
    return out


def summarize_conditions(reports: Dict[str, pd.DataFrame], t_start: float, t_stop: float) -> dict:
    """
    Mean (and min/max for Zd) of the condition time series inside
    [t_start, t_stop] (unix s). Missing reports simply leave their keys out.
    """
    def window(df):
        if df is None or df.empty:
            return None
        sel = df[(df["time"] >= t_start) & (df["time"] <= t_stop)]
        return sel if len(sel) else None

    summary = {}
    cc = window(reports.get("cc"))
    if cc is not None:
        summary.update({
            "cc_source": cc["source"].mode().iloc[0],
            "cc_zd_mean_deg": cc["zd_deg"].mean(),
            "cc_zd_min_deg": cc["zd_deg"].min(),
            "cc_zd_max_deg": cc["zd_deg"].max(),
            "cc_moon_flag": cc["moon_flag"].mode().iloc[0],
        })
    weather = window(reports.get("weather"))
    if weather is not None:
        for col in ("temperature_C", "humidity_pct", "pressure_hPa", "wind_avg_kmh", "wind_peak_kmh"):
            if col in weather:
                summary[col] = weather[col].mean()
        if "wind_peak_kmh" in weather:
            summary["wind_peak_max_kmh"] = weather["wind_peak_kmh"].max()
    pyro = window(reports.get("pyro"))
    if pyro is not None:
        summary["cloudiness_pct"] = pyro["cloudiness_pct"].mean()
        summary["cloudiness_max_pct"] = pyro["cloudiness_pct"].max()
        summary["sky_temperature_K"] = pyro["sky_temperature_K"].mean()
    amc = window(reports.get("amc"))
    if amc is not None:
        idle = amc[amc["amc_panels_2"] == 0]
        if len(idle):
            summary["amc_n_enabled"] = int(idle["amc_panels_0"].median())
            summary["amc_n_undefined"] = int(idle["amc_panels_1"].median())
            summary["amc_n_error"] = int(idle["amc_panels_3"].median())
        summary["amc_n_unfocused"] = amc_unfocused_panels(amc, t_start, t_stop)
    starg = window(reports.get("starg"))
    if starg is not None:
        summary["starg_n_stars"] = starg["starg_n_stars"].mean()
        summary["starg_sky_brightness"] = starg["starg_sky_brightness"].mean()
        summary["starg_dev_rms_deg"] = float(np.sqrt(np.mean(starg["starg_dzd_deg"] ** 2 +
                                                             starg["starg_daz_deg"] ** 2)))
    lidar = window(reports.get("lidar"))
    if lidar is not None:
        summary["lidar_operating"] = bool((lidar["lidar_error_code"] == 0).any())
    return summary


def amc_unfocused_panels(amc: Optional[pd.DataFrame], t_start: float, t_stop: float) -> Optional[int]:
    """
    Number of mirror panels not correctly focused (undefined + error) in
    [t_start, t_stop] (unix s): median over the AMC reports in the window,
    skipping those taken while an adjustment was running. None if there are
    no such reports.
    """
    if amc is None or amc.empty:
        return None
    sel = amc[(amc["time"] >= t_start) & (amc["time"] <= t_stop) & (amc["amc_panels_2"] == 0)]
    if not len(sel):
        return None
    return int(round((sel["amc_panels_1"] + sel["amc_panels_3"]).median()))


def site_location():
    from astropy.coordinates import EarthLocation
    import astropy.units as u
    return EarthLocation(lat=config.ORM_ARRAY.latitude_deg * u.deg,
                         lon=config.ORM_ARRAY.longitude_deg * u.deg,
                         height=config.ORM_ARRAY.height_m * u.m)


def pointing_conditions(ra_deg: float, dec_deg: float, t_start: float, t_stop: float,
                        n_points: int = 20) -> dict:
    """
    Zd/airmass of a target and the moon's altitude, illuminated fraction and
    separation from it, over [t_start, t_stop] (unix s), computed with astropy
    at the ORM site of config.ORM_ARRAY. Available for every night, including
    those without CC reports.
    """
    import astropy.units as u
    from astropy.coordinates import AltAz, SkyCoord, get_body

    times = Time(np.linspace(t_start, t_stop, n_points), format="unix", scale="utc")
    frame = AltAz(obstime=times, location=site_location())
    target = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg).transform_to(frame)
    zd = 90.0 - target.alt.deg
    airmass = 1.0 / np.cos(np.radians(zd))

    moon = get_body("moon", times, site_location())
    sun = get_body("sun", times, site_location())
    elongation = sun.separation(moon).rad
    moon_altaz = moon.transform_to(frame)

    return {
        "zd_mean_deg": float(zd.mean()),
        "zd_min_deg": float(zd.min()),
        "zd_max_deg": float(zd.max()),
        "airmass_mean": float(airmass.mean()),
        "moon_alt_deg": float(moon_altaz.alt.deg.mean()),
        "moon_illumination": float(np.mean((1 - np.cos(elongation)) / 2)),
        "moon_separation_deg": float(moon_altaz.separation(target).deg.mean()),
    }


def conditions_table(log, nights=None) -> pd.DataFrame:
    """
    One row per observation (index night, observation) with its observing
    conditions: Zd/airmass/moon from ``pointing_conditions`` (if the log has
    coordinates), the site-wide CC-report summary (weather, pyrometer,
    starguider, pointing; taken from M2's reports, or the first telescope's),
    the AMC panel counts per telescope as ``<tel>_amc_*``, and the night's /
    observation's ``manual`` overrides. Uses the CSV caches written by
    ``kbo-night prepare``; observations without a .npz cache are skipped.
    """
    rows = []
    for night in log.select(nights):
        cond = load_conditions(night.conditions_dir) if night.conditions_dir else {}
        site_reports = cond.get("m2") or next(iter(cond.values()), {})
        for obs in night.observations.values():
            if not obs.has_cache:
                continue
            t0, t1 = obs.time_span()
            row = {"night": night.name, "observation": obs.name, "type": obs.type}
            if obs.ra_deg is not None:
                row.update(pointing_conditions(obs.ra_deg, obs.dec_deg, t0, t1))
            row.update(summarize_conditions({k: v for k, v in site_reports.items() if k != "amc"}, t0, t1))
            for tel, reports in cond.items():
                amc = summarize_conditions({"amc": reports.get("amc")}, t0, t1)
                row.update({f"{tel}_{k}": v for k, v in amc.items()})
            row.update({k: v for k, v in night.manual_conditions.items() if v is not None})
            row.update(obs.manual)
            rows.append(row)
    return pd.DataFrame(rows).set_index(["night", "observation"])
