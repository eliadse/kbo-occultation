# Functions for data input and processing

import numpy as np
import matplotlib.pyplot as plt
#import astropy.units as u
#from astropy.time import Time

from kbo_occultation import PACKAGE_DATA
from kbo_occultation import config
from kbo_occultation.io import read_stat_binary_file
from kbo_occultation.filtering import highpass_fft


def _channel_signal(data, channel, statistic):
    """
    Extract the flux proxy for one channel from a stat-binary record array. 
    ``"variance"`` (std**2) is proportional to photon flux for shot-noise-dominated 
    AC-coupled PMT signals and is what the multiplicative injection model assumes; 
    ``"std"`` keeps the raw recorded standard deviation.
    """
    std = data[f"std_ch{channel}"]
    if statistic == "variance":
        return std ** 2
    if statistic == "std":
        return std
    raise ValueError(f"statistic must be 'variance' or 'std', got {statistic!r}")


class LightCurve:
    def __init__(self, time, signal, meta=None):
        self.time = time
        self.signal = signal
        self.meta = meta or {}

    @classmethod
    def from_stat_binary(cls,
                     filename,
                     channel,
                     average=None,
                     low_freq_cut=None,
                     time_mode="fixed",
                     sample_time=config.standard_sampling,
                     statistic="variance"):
        """
        Load one channel from a stat binary file.

        - sample_time: The number of DAQ digitizations per statistics record (2**18 in the standard mode); 
        each digitization lasts ``config.standard_sample_duration`` ns, so the record cadence is 
        sample_time * standard_sample_duration ns (~0.5 ms for standard files). 
        - statistic: Selects the flux proxy stored in ``signal`` ("variance" = std**2, "std" = raw std).
        """
        if sample_time is None:
            raise ValueError("A sampling interval (in digitizations) must be provided")

        # Use the binary file function to read an observation file
        data = read_stat_binary_file(filename)

        time = data["time_stamp"].astype(float) / 1.e6 # Now it's in seconds
        t0 = time[0]
        dt = sample_time * config.standard_sample_duration / 1e9
        fs = 1 / dt

        signal = _channel_signal(data, channel, statistic)

        # The time stamps are not reliable, so we re-write them.
        if time_mode == "fixed":
            time = reconstruct_time(len(signal), t0, dt)

        # --- averaging ---
        if average is not None and average > 1:
            signal = average_chunks(signal, average)
            time = average_chunks(time, average)
            fs = fs / average

        # --- filtering ---
        if low_freq_cut is not None:
            signal = highpass_fft(signal, fs, low_freq_cut)

        return cls(time, signal, meta={
            "channel": channel,
            "source": filename,
            "statistic": statistic,
        })

    @classmethod
    def from_stat_binary_all(cls, filename, time_mode="fixed",
                             sample_time=config.standard_sampling,
                             statistic="variance"):
        """
        Load all channels (A, B, C) from a stat binary file as a dict of LightCurves. 
        Same conventions as ``from_stat_binary``; the flux proxy defaults to the variance (std**2).
        """
        if sample_time is None:
            raise ValueError("A sampling interval (in digitizations) must be provided")

        data = read_stat_binary_file(filename)
        time = data["time_stamp"].astype(float) / 1e6
        t0 = time[0]
        dt = sample_time * config.standard_sample_duration / 1e9

        # The time stamps are not reliable, so we re-write them.
        if time_mode == "fixed":
            time = reconstruct_time(len(data["std_chA"]), t0, dt)

        lcs = {}
        for ch in ["A", "B", "C"]:
            lcs[ch] = cls(time, _channel_signal(data, ch, statistic),
                          meta={"channel": ch, "source": filename,
                                "statistic": statistic})
        return lcs
    
    def plot(self, ax=None, **kwargs):
       if ax is None:
           fig, ax = plt.subplots()
       ax.plot(self.time, self.signal, **kwargs)
       return ax

def reconstruct_time(n, t0, dt):
    return t0 + np.arange(n) * dt
    
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