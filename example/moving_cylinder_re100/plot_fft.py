#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.fft import rfft, rfftfreq
from scipy.signal import find_peaks, hilbert
import glob
import argparse
import os

params = {
    'text.usetex': False,
    'font.size': 9,
    'legend.fontsize': 8,
    'legend.handlelength': 2.,
    'agg.path.chunksize': 100000
}
plt.rcParams.update(params)

# Figure output parameters
output_format = 'png'
bbox_option = 'tight'
dpi_quality = 400
fig_width = 3.4
fig_height = 4.3


class LifDrag2D:
    """
    Class to read and process flow data from a file. 
    It detects duplicate lines, reads time series data, and stores them as attributes.
    """
    def __init__(self, filename):
        import numpy as np
        import time

        print("-------------------------------------------------------")
        print("Reading " + filename)
        start_time = time.perf_counter()

        last_line = None
        data_lines = []
        with open(filename, "r") as f:
            for line in f:
                parts = line.strip().split()
                line_without_first = " ".join(parts[1:])
                if line_without_first != last_line:
                    data_lines.append(list(map(float, parts)))
                else:
                    print(f"Skipped duplicate line: {line.strip()}")
                last_line = line_without_first

        data_array = np.transpose(data_lines)

        expected_attributes = ["i", "t", "dgx", "dpx", "dvx", "dgy", "dpy", "dvy", "tqz", "tpz", "tvz"]
        if data_array.shape[0] != len(expected_attributes):
            print(f"Warning: Number of columns ({data_array.shape[0]}) differs from expected ({len(expected_attributes)}).")
            print("Using generic attribute names.")
            attributes = [f"col_{i}" for i in range(data_array.shape[0])]
        else:
            attributes = expected_attributes

        for idx, attr in enumerate(attributes):
            setattr(self, attr, data_array[idx])

        del data_array, data_lines
        end_time = time.perf_counter()
        print(f"Reading and processing took {end_time - start_time:.4f} s.")


def compute_fft(signal, dt, normalize=False, print_peaks=True):
    """
    Compute FFT and PSD of a given signal, detect top peaks, and optionally print them.
    """
    import time
    start_time = time.perf_counter()

    N = len(signal)
    fft_vals = rfft(signal)
    freqs = rfftfreq(N, dt)
    psd = (np.abs(fft_vals)**2) * (dt / N)

    if normalize:
        # Adjust normalization if needed
        psd /= (N / dt)

    # Detect peaks
    threshold = np.min(psd)
    peak_indices = find_peaks(psd, height=threshold, prominence=threshold)[0]

    num_to_print = 4
    top_n = min(num_to_print, len(peak_indices))
    if top_n > 0:
        partial_indices = np.argpartition(-psd[peak_indices], top_n - 1)[:top_n]
        sorted_partial = partial_indices[np.argsort(-psd[peak_indices][partial_indices])]
        top_freqs = freqs[peak_indices][sorted_partial]
        top_psd_values = psd[peak_indices][sorted_partial]

        if print_peaks:
            for i in range(top_n):
                f_val = top_freqs[i]
                omega_val = 2 * np.pi * f_val
                print("PSD: {:.4f}, f: {:.8f}, ω: {:.8f}".format(top_psd_values[i], f_val, omega_val))
    else:
        print("No peaks found.")

    end_time = time.perf_counter()
    print(f"FFT took {end_time - start_time:.4f} s.")

    prfreq = top_freqs if top_n > 0 else np.array([])
    return psd, freqs, prfreq, (top_psd_values if top_n > 0 else np.array([]))


def plot_probe_data():
    # Input files
    files = glob.glob('*.his')

    # Set time trimming parameters
    time_skips = [0,0]  # set to 0 to use the full time range
    time_maxes = [0,0]  # set to 0 to use the full time range

    for i, filename in enumerate(files):
        print(f'Opening file {filename}')

        # Get root name from input file
        root_name = filename.split('.')[0]

        # Read number of probes
        with open(filename, 'r') as file:
            n_probes = int(file.readline().strip())
            print(f'Number of probes found: {n_probes}')

            # Read probe coordinates
            coords = np.zeros((n_probes, 3))
            for n in range(n_probes):
                line = file.readline().strip()
                coords[n] = line.split()

        # Load the remaining data
        data = np.loadtxt(filename, skiprows=n_probes + 1)
        total_data_lines, num_columns = data.shape

        if total_data_lines % n_probes != 0:
            raise ValueError("Data lines are not divisible by the number of probes. File format may be inconsistent.")

        # Reshape data into (time_steps, n_probes, num_columns)
        time_steps = total_data_lines // n_probes
        data = data.reshape(time_steps, n_probes, num_columns)

        print(f"Data shape: {data.shape}")
        print(f"Coordinates shape: {coords.shape}")
        print(f"Columns per probe: {num_columns}")
        print(f"Number of time steps: {time_steps}")

        for probe_idx in range(n_probes):
            print(f'Probe number {probe_idx + 1}')

            # Extract time and velocities from data
            t = data[:, probe_idx, 0][1:]
            v_comp = data[:, probe_idx, 2][1:]  # vertical velocity

            # Determine trimming times
            t_min_raw, t_max_raw = t.min(), t.max()
            print(f'Original time series: {t_min_raw} to {t_max_raw}')

            t_min = max(t_min_raw, float(time_skips[i])) if float(time_skips[i]) > 0 else t_min_raw
            t_max = min(t_max_raw, float(time_maxes[i])) if float(time_maxes[i]) > 0 else t_max_raw
            print(f'Trimmed time series: {t_min} to {t_max}')

            itmin, itmax = np.where(t >= t_min)[0][0], np.where(t >= t_max)[0][0]
            t_crop = t[itmin:itmax]
            dt_array = np.zeros(len(t_crop), dtype=np.float64)
            for xx in range(1, len(dt_array) - 1):
                dt_array[xx] = t_crop[xx + 1] - t_crop[xx]

            dt = dt_array.max()
            t_new = np.linspace(t[itmin], t[itmax], int((t[itmax] - t[itmin]) / dt), endpoint=False)

            fig, axs = plt.subplots(2, sharex=False, figsize=(fig_height, fig_width))
            axs[0].set_xlabel(r'$t$')
            axs[0].set_title(f'x,y,z = {coords[probe_idx,0]:.2f}, {coords[probe_idx,1]:.2f}, {coords[probe_idx,2]:.2f}')

            # Process vertical velocity
            v_mean = np.mean(v_comp[itmin:itmax])
            v_adj = v_comp[itmin:itmax] - v_mean
            print('(v) mean, min, max =', round(v_mean, 6), v_adj.min(), v_adj.max())

            v_interp = interpolate.interp1d(t_crop, v_adj, kind='slinear', fill_value="extrapolate")(t_new)

            # Reference frequency for plotting
            ref_freq = 0.1643
            period = 1.0 / ref_freq
            period2 = 1.0 / (ref_freq / 2)
            axs[1].axvline(ref_freq, ls='dashed', lw=0.6, c='k', label='$St=$'+("{0:.4f}".format(round(ref_freq,4))))
            axs[1].axvline(ref_freq/2, ls='dashed', lw=0.6, c='blue', label='$St/2$')

            # Plot the time series
            axs[0].axhline(0.00, ls='--', lw=0.5, c='k')
            axs[0].scatter(t_new, v_interp, c='r', s=0.1)

            axs[0].axvline(period, ls='--', lw=0.5)
            n_cycles = int(max(t_new) / period)
            for cycle_idx in range(1, n_cycles + 1):
                axs[0].axvline(cycle_idx * period, ls='--', lw=0.1, c='gray')

            axs[0].set_ylabel(r'$v^{\prime}$')

            # Compute and plot FFT
            psd, freqs, peak_freqs, peak_psd_vals = compute_fft(v_interp, dt)
            axs[1].plot(freqs, np.abs(psd), c='r', lw=0.6)
            if len(peak_freqs) > 0:
                axs[1].axvline(x=peak_freqs[0], c='r', lw=0.5, ls='--', label='$St=$'+("{0:.4f}".format(round(peak_freqs[0],4))))
                print('Number of cycles: {:.2f}'.format((t_new[-1] - t_new[0]) / (1 / peak_freqs[0])))

            axs[1].set_xlabel(r'$St$')
            axs[1].set_xscale('log')
            axs[1].set_yscale('log')

            axs[1].set_xlim(5e-2, 4.0)
            axs[1].set_ylim(bottom=1e-7)

            plt.legend(loc='upper right')

            fname = f'{root_name}_his{probe_idx+1}_fft.{output_format}'
            print('Saving', fname)
            plt.savefig(fname, format=output_format, dpi=dpi_quality, bbox_inches=bbox_option)
            print('------------------------------------------')
            plt.close()


def plot_drag_data():
    # Input file and parameters
    filename = './dragxy.dat'
    data = LifDrag2D(filename)

    # Time series parameters
    t_min = 0.1  # Starting time
    t_max = 0  # Ending time (0 means use all data)
    skip = 2
    t = data.t[skip:]  # Time array, skip first point
    
    # Skip first point for all data
    dgx = data.dgx[skip:]
    dgy = data.dgy[skip:]

    # Determine time range
    if t_max <= 0:
        t_max = t[-1]
    print(f'Time series: {t_min} to {t_max}')

    # Find indices for time range
    itmin = np.where(t >= t_min)[0][0]
    itmax = np.where(t >= t_max)[0][0] if t_max < t[-1] else len(t) - 1
    t_crop = t[itmin:itmax]

    # Calculate time step
    dt_array = np.zeros(len(t_crop), dtype=np.float64)
    for i in range(1, len(dt_array) - 1):
        dt_array[i] = t_crop[i + 1] - t_crop[i]
    dt = dt_array.max()

    # Create new time array
    t_new = np.linspace(t[itmin], t[itmax], int((t[itmax] - t[itmin]) / dt), endpoint=False)

    # Process drag forces
    dgx_crop = dgx[itmin:itmax]
    dgy_crop = dgy[itmin:itmax]

    # Interpolate drag forces
    dgx_interp = interpolate.interp1d(t_crop, dgx_crop, kind='slinear', fill_value="extrapolate")(t_new)
    dgy_interp = interpolate.interp1d(t_crop, dgy_crop, kind='slinear', fill_value="extrapolate")(t_new)

    # Create figure
    fig, axs = plt.subplots(2, sharex=False, figsize=(fig_height, fig_width))

    # Plot time series
    axs[0].set_xlabel(r'$t$')
    axs[0].set_ylabel(r'$F_x,F_y$')
    axs[0].plot(t_new, dgx_interp, c='r', label=r'$F_x$')
    axs[0].plot(t_new, dgy_interp, c='b', label=r'$F_y$')
    axs[0].legend(loc='upper right')

    # Reference frequency for plotting
    ref_freq = 0.1643
    period = 1.0 / ref_freq
    period2 = 1.0 / (ref_freq / 2)

    # Plot period markers
    axs[0].axvline(period, ls='--', lw=0.5)
    n_cycles = int(max(t_new) / period)
    for cycle_idx in range(1, n_cycles + 1):
        axs[0].axvline(cycle_idx * period, ls='--', lw=0.1, c='gray')

    # Compute and plot FFT for both components
    psd_x, freqs_x, peak_freqs_x, _ = compute_fft(dgx_interp, dt)
    psd_y, freqs_y, peak_freqs_y, _ = compute_fft(dgy_interp, dt)

    axs[1].plot(freqs_x, np.abs(psd_x), c='r', lw=0.6, label=r'$F_x$')
    axs[1].plot(freqs_y, np.abs(psd_y), c='b', lw=0.6, label=r'$F_y$')

    # Plot reference frequencies
    axs[1].axvline(ref_freq, ls='dashed', lw=0.6, c='k', label='$St=$'+("{0:.4f}".format(round(ref_freq,4))))
    axs[1].axvline(ref_freq/2, ls='dashed', lw=0.6, c='gray', label='$St/2$')

    if len(peak_freqs_x) > 0:
        axs[1].axvline(x=peak_freqs_x[0], c='r', lw=0.5, ls='--')
    if len(peak_freqs_y) > 0:
        axs[1].axvline(x=peak_freqs_y[0], c='b', lw=0.5, ls='--')

    axs[1].set_xlabel(r'$St$')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlim(5e-2, 4.0)
    axs[1].set_ylim(bottom=1e-7)
    axs[1].legend(loc='upper right')

    # Save figure
    fname = 'drag_fft.png'
    print('Saving', fname)
    print('------------------------------------------')

    plt.savefig(fname, format=output_format, dpi=dpi_quality, bbox_inches=bbox_option)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot FFT of probe or drag data')
    parser.add_argument('--drag', action='store_true', help='Plot drag data instead of probe data')
    args = parser.parse_args()

    # Always try to plot drag data first if file exists
    #drag_file = './dragxy.dat'
    #if os.path.exists(drag_file):
    #    plot_drag_data()
    
    # Then plot probe data unless --drag was specified
    if not args.drag:
        plot_probe_data()
