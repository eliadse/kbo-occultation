# Functions for data input and processing

import numpy as np
import matplotlib.pyplot as plt
#import astropy.units as u
#from astropy.time import Time

from kbo_occultation import PACKAGE_DATA
from kbo_occultation.io import read_stat_binary_file

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
                     fs=None,
                     time_mode="fixed",
                     sample_time=262144):
       
        # Use the binary file function to read an observation file
        data = read_stat_binary_file(filename)

        # In the standard, the statistics of 2^18 digitizations are calculated
        # This corresponds to ~0.5ms (the DAQ digitizes at 500MSamples/s --> 2ns)
        # The fast one takes 2^17 digitizations, corresponding to ~0.1ms
        #TODO Add the standard sampling and standard sample size to a config or something
        # so that it's easy to know what I'm using. And use it from the config, not hard coded
        sample_time = sample_time * 2 # The 2 is for the 2ns

        #time = Time(np.array(data['time_stamp']) / 1.e6, format='unix', scale='utc')
        time = data["time_stamp"].astype(float) / 1.e6 # Now it's in seconds
        t0 = time[0]
        print("t0 = ", t_0)
        #time -= t0
        dt = sample_time / 1e9
        fs = 1/dt
        
        signal = data[f"std_ch{channel}"]**2
        
        # The time stamps are not reliable, so we re-write them.
        if time_mode == "fixed":
            if sample_time is None:
                raise ValueError("A sampling interval must be provided when fixed_time=True")

            time = reconstruct_time(len(signal), t0, dt)
        
        # --- averaging ---
        if average is not None and average > 1:
            signal = average_chunks(signal, average)
            time = average_chunks(time, average)

        # --- filtering ---
        if low_freq_cut is not None:
            signal = low_frequency_cut(signal, fs, low_freq_cut)

        return cls(time, signal, meta={
            "channel": channel,
            "source": filename
        })

    @classmethod
    def from_stat_binary_all(cls, filename, time_mode="fixed", sample_time=262144):
        data = read_stat_binary_file(filename)
        time = data["time_stamp"].astype(float) / 1e6
        t0 = time[0]

        sample_time = sample_time * 2 # The 2 is for the 2ns

        # The time stamps are not reliable, so we re-write them.
        if time_mode == "fixed":
            if sample_time is None:
                raise ValueError("A sampling interval must be provided when fixed_time=True")

            time = reconstruct_time(len(data["std_chA"]), t0, sample_time / 1e9)

        lcs = {}
        for ch in ["A", "B", "C", "D"]:
            lcs[ch] = cls(time, data[f"std_ch{ch}"],
                          meta={"channel": ch, "source": filename})
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

    def average_chunks(x, y, n_samples):
        # Take the array, remove the last points which go over lenx/n_samples
        # Reshape the array so that rows contain n_samples columns.
        # Take the average of the rows.
        averaged_x = np.average(x[0:-(len(x)%n_samples)].reshape(-1, n_samples), axis=1)
        averaged_y = np.average(y[0:-(len(y)%n_samples)].reshape(-1, n_samples), axis=1)
        error_y = np.std(y[0:-(len(y)%n_samples)].reshape(-1, n_samples), axis=1)/np.sqrt(n_samples)
        return averaged_x, averaged_y, error_y

#minimalistic
def average_chunks(x, n):
    if n <= 1:
        return x
    return np.mean(x[:len(x)//n*n].reshape(-1, n), axis=1)


def low_frequency_cut(y_values, frequency_cut):
        from scipy.fftpack import fftfreq, irfft, rfft
        dt = sample_time.to('s').value
        w = fftfreq(y_values.size, d=dt)
        f_signal = rfft(y_values)
        cut_f_signal = f_signal.copy()
        # cut signal below frequency_cut
        cut_f_signal[(np.abs(w) < frequency_cut)] = 0
        cs = irfft(cut_f_signal)
        return cs

#minimalistic
def low_frequency_cut(y, fs, fcut):
    freqs = np.fft.rfftfreq(len(y), d=1/fs)
    fft = np.fft.rfft(y)

    fft[np.abs(freqs) < fcut] = 0
    return np.fft.irfft(fft, n=len(y))

#def load_star_info()

def plot_test_stars(data_dict):
    fig = plt.figure(figsize=[12,7])
    #for star in data_dict:
    for star in mag_dict:
        data = data_dict[star][0]
        # mask for outliers
        #mask = data < (np.mean(data) + 4* np.std(data))
        # For tests
        #plt.scatter(data_dict[star][1], data, marker=".", label=f'{star}', alpha=0.7)
        plt.scatter(data_dict[star][1], data, marker=".", label=f'{star}, Mag(B): {mag_dict[star]}', alpha=0.7)
        print(star, np.mean(data))

    plt.legend(loc="lower right")
    plt.xlabel("Time [s]")
    plt.ylabel(r"std_dev$^2$")
    plt.title('Fast photometry of test stars with MAGIC-II')
    plt.tight_layout()
    #plt.savefig(f'{data_path}/All_together.png')
    plt.show()