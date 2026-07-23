# Functions for data input and processing

import json
import os

import numpy as np
import matplotlib.pyplot as plt
#import astropy.units as u
#from astropy.time import Time

from kbo_occultation import PACKAGE_DATA
from kbo_occultation import config
from kbo_occultation.io import read_stat_binary_file
from kbo_occultation.filtering import highpass_fft


def _channel_signal(data, channel):
    """
    Extract the flux proxy for one channel from a stat-binary record array.
    Always the variance (std**2): it is proportional to photon flux for
    shot-noise-dominated AC-coupled PMT signals and is what the multiplicative
    injection model assumes.
    """
    return data[f"std_ch{channel}"] ** 2


class LightCurve:
    def __init__(self, time, signal, meta=None):
        self.time = time
        self.signal = signal
        self.meta = meta or {}
        # The recorded timestamps are unreliable, so the loaders always rebuild
        # the time grid at a fixed cadence (reconstruct_time). dt is therefore a
        # constant read straight off the grid, not something to re-derive from
        # the timestamps at each use.
        self.dt = float(time[1] - time[0])

    @classmethod
    def from_stat_binary(cls,
                     filename,
                     channel,
                     average=None,
                     low_freq_cut=None,
                     sample_time=config.standard_sampling):
        """
        Load one channel from a stat binary file.

        - sample_time: The number of DAQ digitizations per statistics record (2**18 in the standard mode);
        each digitization lasts ``config.standard_sample_duration`` ns, so the record cadence is
        sample_time * standard_sample_duration ns (~0.5 ms for standard files).

        The recorded timestamps are unreliable, so the time grid is always
        rebuilt at this fixed cadence; the flux proxy stored in ``signal`` is
        always the variance (std**2).
        """
        #TODO note that from_stat_binary and from_stat_binary_all are different. You should check if
        # this affects something down the line
        if sample_time is None:
            raise ValueError("A sampling interval (in digitizations) must be provided")

        # Use the binary file function to read an observation file
        data = read_stat_binary_file(filename)

        t0 = data["time_stamp"][0].astype(float) / 1.e6  # Now it's in seconds
        dt = sample_time * config.standard_sample_duration / 1e9
        fs = 1 / dt

        signal = _channel_signal(data, channel)

        # The time stamps are not reliable, so we always re-write them.
        time = reconstruct_time(len(signal), t0, dt)

        time, signal = _finish_signal(time, signal, fs, average, low_freq_cut)

        return cls(time, signal, meta={
            "channel": channel,
            "source": filename,
        })

    @classmethod
    def from_stat_binary_all(cls, filename,
                             sample_time=config.standard_sampling):
        """
        Load all channels (A, B, C) from a stat binary file as a dict of
        LightCurves. Same conventions as ``from_stat_binary``: the time grid is
        always rebuilt at a fixed cadence and the flux proxy is the variance
        (std**2).
        """
        if sample_time is None:
            raise ValueError("A sampling interval (in digitizations) must be provided")

        data = read_stat_binary_file(filename)
        t0 = data["time_stamp"][0].astype(float) / 1e6
        dt = sample_time * config.standard_sample_duration / 1e9

        # The time stamps are not reliable, so we always re-write them.
        time = reconstruct_time(len(data["std_chA"]), t0, dt)

        lcs = {}
        for ch in ["A", "B", "C"]:
            lcs[ch] = cls(time, _channel_signal(data, ch),
                          meta={"channel": ch, "source": filename})
        return lcs

    @classmethod
    def from_preprocessed(cls, path, channel, *, cleaned=True,
                          average=None, low_freq_cut=None):
        """
        Load one channel from a pre-processed ``.npz`` cache written by
        ``preprocess_observation``.

        The cache stores, per channel, both the raw variance (``var_{ch}_raw``)
        and the outlier-cleaned variance (``var_{ch}_clean``); ``cleaned``
        selects between them (default: the cleaned series). The uniform time
        grid is rebuilt from the stored ``t0``/``dt``/``n`` with the same
        ``reconstruct_time`` convention as the binary loaders, and the optional
        ``average``/``low_freq_cut`` post-processing is identical to
        ``from_stat_binary``.
        """
        with np.load(path, allow_pickle=False) as npz:
            meta = json.loads(str(npz["meta"]))
            t0 = float(npz["t0"])
            dt = float(npz["dt"])
            n = int(npz["n"])
            key = f"var_{channel}_clean" if cleaned else f"var_{channel}_raw"
            signal = np.asarray(npz[key], dtype=float)

        fs = 1 / dt
        time = reconstruct_time(n, t0, dt)
        time, signal = _finish_signal(time, signal, fs, average, low_freq_cut)

        return cls(time, signal, meta={
            "channel": channel,
            "source": meta.get("source", path),
            "cleaned": cleaned,
            "outlier_threshold_sigma": meta.get("outlier_threshold_sigma"),
            "outlier_window_samples": meta.get("outlier_window_samples"),
            "outlier_polarity": meta.get("outlier_polarity"),
        })

    @classmethod
    def from_preprocessed_all(cls, path, *, cleaned=False):
        """
        Load all channels (A, B, C) from a pre-processed ``.npz`` cache as a
        dict of LightCurves -- the analogue of ``from_stat_binary_all``.

        ``cleaned`` defaults to ``False`` so that pointing an existing raw-data
        workflow (e.g. the search pipeline) at a cache reproduces the previous
        raw-variance behaviour exactly; callers opt into the outlier-cleaned
        series explicitly.
        """
        return {ch: cls.from_preprocessed(path, ch, cleaned=cleaned)
                for ch in ["A", "B", "C"]}

    def plot(self, ax=None, **kwargs):
       if ax is None:
           fig, ax = plt.subplots()
       ax.plot(self.time, self.signal, **kwargs)
       return ax

def reconstruct_time(n, t0, dt):
    return t0 + np.arange(n) * dt


def _finish_signal(time, signal, fs, average, low_freq_cut):
    """
    Apply the optional block-averaging and high-pass steps shared by the
    binary and pre-processed loaders. Returns the (possibly shortened)
    ``(time, signal)`` pair.
    """
    # TODO check if the order of averaging and filtering is relevant
    # --- averaging ---
    if average is not None and average > 1:
        signal = average_chunks(signal, average)
        time = average_chunks(time, average)
        fs = fs / average

    # --- filtering ---
    if low_freq_cut is not None:
        signal = highpass_fft(signal, fs, low_freq_cut)

    return time, signal


def preprocess_observation(bin_path, out_path=None, *,
                           sample_time=config.standard_sampling,
                           outlier_kwargs=None):
    """
    Pre-process a raw stat-binary observation into a compact ``.npz`` cache.

    Reads the binary once, extracts only the information the package actually
    uses (the per-channel variance flux proxy and the start time), applies the
    outlier correction, and stores both the raw and the cleaned variance per
    channel so that nothing downstream changes silently. The reconstructed time
    grid is captured as ``t0``/``dt``/``n`` rather than a full array.

    Parameters
    ----------
    bin_path : str
        Path to the raw stat-binary observation file.
    out_path : str, optional
        Destination ``.npz`` path. Defaults to ``bin_path`` with its extension
        replaced by ``.npz``.
    sample_time : int
        Digitizations per statistics record; sets the record cadence, exactly
        as in ``from_stat_binary`` (default: ``config.standard_sampling``).
    outlier_kwargs : dict, optional
        Extra keyword arguments forwarded to ``remove_outliers_lightcurve``
        (e.g. ``threshold_sigma``, ``window_samples``, ``polarity``).

    Returns
    -------
    str
        The path of the written ``.npz`` file.
    """
    from kbo_occultation.filtering import remove_outliers_lightcurve

    if sample_time is None:
        raise ValueError("A sampling interval (in digitizations) must be provided")

    outlier_kwargs = dict(outlier_kwargs or {})

    data = read_stat_binary_file(bin_path)
    t0 = data["time_stamp"][0].astype(float) / 1e6
    dt = sample_time * config.standard_sample_duration / 1e9
    n = len(data["std_chA"])
    time = reconstruct_time(n, t0, dt)

    arrays = {"t0": np.float64(t0), "dt": np.float64(dt), "n": np.int64(n)}
    outlier_meta = {}
    for ch in ["A", "B", "C"]:
        var_raw = _channel_signal(data, ch)
        cleaned_lc = remove_outliers_lightcurve(
            LightCurve(time, var_raw), **outlier_kwargs)
        arrays[f"var_{ch}_raw"] = var_raw.astype(np.float32)
        arrays[f"var_{ch}_clean"] = cleaned_lc.signal.astype(np.float32)
        # The outlier parameters are identical across channels; record once.
        outlier_meta = {
            "outlier_threshold_sigma": cleaned_lc.meta.get("outlier_threshold_sigma"),
            "outlier_window_samples": cleaned_lc.meta.get("outlier_window_samples"),
            "outlier_polarity": cleaned_lc.meta.get("outlier_polarity"),
        }

    meta = {
        "source": os.path.basename(bin_path),
        "sample_time": int(sample_time),
        "dt": float(dt),
        "channels": ["A", "B", "C"],
        **outlier_meta,
    }
    arrays["meta"] = json.dumps(meta)

    if out_path is None:
        out_path = os.path.splitext(bin_path)[0] + ".npz"
    np.savez_compressed(out_path, **arrays)
    return out_path


def plot_lightcurves(lightcurves, labels=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    for i, lc in enumerate(lightcurves):
        label = None
        if labels is not None:
            label = labels[i]
        elif "channel" in lc.meta:
            label = f"ch{lc.meta['channel']}"

        ax.plot(lc.time, lc.signal, label=label)

    ax.set_xlabel("Time")
    ax.set_ylabel("Signal")

    if labels is not None or any("channel" in lc.meta for lc in lightcurves):
        ax.legend()

    return ax


# TODO Add a function to set the name and magnitude of the star. Maybe its size too

def average_chunks(x, n):
    """
    Trim x to a multiple of n and average non-overlapping blocks of n.
    See also detectability.bin_average, which does the same block
    averaging on a (time, signal) pair and propagates sigma.
    """
    if n <= 1:
        return x
    return np.mean(x[:len(x)//n*n].reshape(-1, n), axis=1)