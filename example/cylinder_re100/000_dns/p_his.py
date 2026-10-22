#!/usr/bin/env python
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
import os, os.path, csv, time, pickle, glob

from scipy import signal
from scipy.interpolate import interp1d

try:
    from scipy.fft import rfft, rfftfreq  # SciPy >= 1.4.0
except ImportError:
    from scipy.fftpack import rfft, rfftfreq  # SciPy < 1.4.0


def periodogram_rfft(x, fs):
    """Compute PSD using periodogram with real FFT."""
    freqs, psd = signal.periodogram(x, fs, scaling="spectrum")
    return freqs, psd


def find_peaks(st, psd, threshold=0.01):
    """Return peak locations and values above a fraction of the maximum."""
    if len(psd) == 0:
        return np.array([]), np.array([])
    peak_indices = signal.find_peaks(psd, height=max(psd) * threshold)[0]
    return st[peak_indices], psd[peak_indices]


params = {"text.usetex": False, "font.size": 11, "legend.fontsize": 11, "legend.handlelength": 2.5, "agg.path.chunksize": 100000}
plt.rcParams.update(params)

formt = "png"
ajust = "tight"
qual = 400
fig_width = 4.3
fig_height = 16 * fig_width / 9

files = [
    "1cyl.his",
]
tskps = [0]  # set to 0 to plot all
tmaxs = [0]

for i, filen in enumerate(files):
    print(f"Opening file {filen}")

    with open(filen, "r") as file:
        nps = int(file.readline().rstrip())
        print(f"Number of probes found: {nps}")

        coords = np.zeros((nps, 3))
        for n, line in zip(range(nps), file):
            coords[n] = line.split()

    # careful with race condition if the simulation is not finished
    # we might "accidentally" read negative tmin.
    max_retries = 3
    retry_count = 0
    data_valid = False
    
    while not data_valid and retry_count < max_retries:
        try:
            raw_data = np.loadtxt(filen, skiprows=nps + 1)
            
            # Try to determine if it's 2D or 3D by attempting to reshape
            try:
                data = raw_data.reshape(-1, nps, 5)
                is_3d = True
            except:
                data = raw_data.reshape(-1, nps, 4)
                is_3d = False
                
            # Verify data integrity by checking the first probe's time values
            t = data[:, 0, 0][1:]  # Check first probe's time values
            if t.size > 0 and t.min() >= 0:
                data_valid = True
            else:
                retry_count += 1
                if retry_count < max_retries:
                    time.sleep(0.1)  # Small delay before retry
                    
        except (ValueError, IndexError):
            retry_count += 1
            if retry_count < max_retries:
                time.sleep(0.1)  # Small delay before retry
    
    if not data_valid:
        print(f"Warning: Could not read valid data from {filen} after {max_retries} attempts. Skipping.")
        continue

    for prb in range(0, nps, 1):
        print(f"Probe number {prb + 1}")

        t = data[:, prb, 0][1:]
        u = data[:, prb, 1][1:]
        v = data[:, prb, 2][1:]
        if is_3d:
            w = data[:, prb, 3][1:]
            p = data[:, prb, 4][1:]
        else:
            p = data[:, prb, 3][1:]

        # Adjust time series (crop and interpolate for constant time step)
        tmin, tmax = t.min(), t.max()
        print(f"Original time series: {tmin} {tmax}")

        tmin = max(tmin, float(tskps[i])) if float(tskps[i]) > 0 else tmin
        tmax = min(tmax, float(tmaxs[i])) if float(tmaxs[i]) > 0 else tmax
        print(f"Cropped time series: {tmin} {tmax}")

        itmin, itmax = np.where(t >= tmin)[0][0], np.where(t >= tmax)[0][0]
        to = t[itmin:itmax]
        dtt = np.zeros(len(to), dtype=np.float64)
        for xx in range(1, len(dtt) - 1):
            dtt[xx] = to[xx + 1] - to[xx]

        dt = dtt.max()
        tn = np.linspace(t[itmin], t[itmax], int((t[itmax] - t[itmin]) / dt), endpoint=False)

        # Figure
        fig, axs = plt.subplots(2, sharex=False)
        fig.set_size_inches(fig_height, fig_width)

        # Plot time series (signal on upper part of the figure)
        axs[0].set_xlabel(r"$t$")
        axs[0].set_ylabel(r"$v$")
        axs[0].set_title(f"x,y,z={coords[prb][0]},{coords[prb][1]},{coords[prb][2]}")

        # We focus on the vertical velocity signal
        nm = "v"  # name of the signal
        cor = "r"  # color of the signal
        vo = v[itmin:itmax]
        qur = np.mean(vo)  # compute mean
        vo -= qur  # remove mean
        print("(", nm, ") mean, min, max=", round(qur, 6), vo.min(), vo.max())
        vn = interp1d(to, vo, kind="nearest", fill_value="extrapolate")(tn)  # snap to nearest value

        axs[0].scatter(tn, vn, c=cor, ls="-", s=0.001)  # plot signal

        # Compute PSD using periodogram_rfft
        ufreq, psd = periodogram_rfft(vn, fs=1.0 / dt)

        # Find peaks in the PSD
        peak_freqs, peak_vals = find_peaks(ufreq, psd, threshold=0.01)

        # Plot PSD
        axs[1].semilogy(ufreq, psd, c=cor, lw=0.6)
        axs[1].scatter(peak_freqs, peak_vals, c="k", marker="x", s=40, label="peaks")
        if len(peak_freqs) > 0:
            axs[1].axvline(x=peak_freqs[0], c=cor, lw=0.5, ls="--", label=f"$St={peak_freqs[0]:.4f}$")
        axs[1].set_xlabel(r"$St$")
        axs[1].set_xlim(1e-4, 1)
        axs[1].axvline(x=0.125, c="k", lw=0.5, ls="--", label=f"Ref. $Re=50$ $St=0.125$")

        psd_without_0 = psd[4:]
        axs[1].set_ylim(psd_without_0.min(), psd_without_0.max() * 10)
        plt.legend(loc="upper right")
        file_root = filen.split(".")[0]
        fname = file_root + str(prb + 1) + "_fft." + formt
        print("Saving ", fname)
        print()
        plt.savefig(fname, format=formt, dpi=qual, bbox_inches=ajust)
        plt.close()
        plt.clf()

        # Phase space plot
        fig = plt.figure(figsize=(fig_height, fig_width))
        plt.xlabel(r"$v$")
        plt.ylabel(r"$\dot{v}$")
        plt.plot(v[2:-2], np.gradient(v)[2:-2] / dt, c="r", ls="-", lw=0.8)

        fname = file_root + str(prb + 1) + "_phase_space." + formt
        print("Saving ", fname)
        plt.savefig(fname, format=formt, dpi=qual, bbox_inches=ajust)
