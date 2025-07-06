#!/usr/bin/env python
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os, sys, argparse, shutil

try:
    sys.path.insert(0, os.path.join(os.environ["NEKSTAB_SOURCE_ROOT"], "bin"))
except KeyError:
    raise EnvironmentError("NEKSTAB_SOURCE_ROOT environment variable is not set.")
from nekStab_tools import (
    interpolate_signal,
    periodogram_rfft,
    find_peaks,
    LiftDragLoader,
    get_case_name_from_usr,
)

# --- Plotting Parameters ---
plot_params = {
    "text.usetex": shutil.which("latex") is not None,
    "font.size": 11,
    "legend.fontsize": 11,
    "legend.handlelength": 2.5,
    "agg.path.chunksize": 100000,
    "format": "png",
    "adjust": "tight",
    "quality": 400,
    "fig_width": 4.3,
    "fig_height": 16 * 4.3 / 9,
}

# Update rcParams with only the valid keys to avoid errors
valid_rc_keys = [
    "text.usetex",
    "font.size",
    "legend.fontsize",
    "legend.handlelength",
    "agg.path.chunksize",
]
rc_params_to_set = {key: plot_params[key] for key in valid_rc_keys}
plt.rcParams.update(rc_params_to_set)


# --- Plotting Functions ---
def plot_fft(ax, tn, vn, dt, label, color, scaling="spectrum", plot_reference=False):
    """Computes and plots the FFT (periodogram) of a signal."""
    ufreq, psd = periodogram_rfft(vn, fs=1.0 / dt, scaling=scaling)
    peak_freqs, peak_vals = find_peaks(ufreq, psd, threshold=0.01)

    # Plot reference line first so it appears first in legend
    if plot_reference:
        ax.axvline(x=0.125, c="k", lw=0.5, ls="--", label=f"Ref. $Re=50$ $St=0.125$")

    ax.loglog(ufreq, psd, c=color, lw=0.8, label=label)
    ax.scatter(peak_freqs, peak_vals, c="k", marker="x", s=40)
    if len(peak_freqs) > 0:
        # Find the highest peak instead of the first peak
        highest_peak_idx = np.argmax(peak_vals)
        highest_peak_freq = peak_freqs[highest_peak_idx]
        ax.axvline(x=highest_peak_freq, c=color, lw=0.5, ls="--", label=f"$St={highest_peak_freq:.4f}$")

    nyquist = 0.5 / dt
    ax.set_xlim(1e-2, nyquist)
    psd_without_0 = psd[4:]
    if len(psd_without_0) > 0:
        ax.set_ylim(psd_without_0.min(), psd_without_0.max() * 10)
    ax.set_xlabel(r"$St$")
    if scaling == "spectrum":
        ax.set_ylabel("Power Spectrum")
    elif scaling == "density":
        ax.set_ylabel("PSD")
    else:
        ax.set_ylabel("Power")
    # Legend and title are handled externally


def plot_phase_portrait(ax, vn, dt, color, name="v", add_label=False):
    """Computes and plots the phase portrait of a signal."""
    vdot_interp = np.gradient(vn, dt)
    vdot_trim = vdot_interp[2:-2]
    v_trim = vn[2:-2]

    ax.set_xlabel(f"${name}$")
    ax.set_ylabel(rf"$\dot{{{name}}}$")
    ax.set_title("Phase Portrait (trimmed)")
    if add_label:
        ax.plot(v_trim, vdot_trim, c=color, lw=0.8, label=f"${name}$")
    else:
        ax.plot(v_trim, vdot_trim, c=color, lw=0.8)


# --- Data Processors ---
def process_his_file(his_file):
    """Load, process, and plot data from a .his file."""
    print(f"--- Processing History File: {his_file} ---")
    with open(his_file, "r") as file:
        try:
            nps = int(file.readline().rstrip())
        except (ValueError, IndexError):
            print(f"Error: Could not read number of probes from {his_file}. Skipping.")
            return
        print(f"Number of probes found: {nps}")
        coords = np.zeros((nps, 3))
        for n, line in zip(range(nps), file):
            try:
                coord_parts = line.split()
                if len(coord_parts) >= 3:
                    coords[n] = coord_parts[:3]  # Take first 3 coordinates
                else:
                    coords[n] = [0, 0, 0]  # Default coordinates if missing
            except (ValueError, IndexError):
                coords[n] = [0, 0, 0]  # Default coordinates on error

    data_lines = []
    with open(his_file, "r") as file:
        lines = file.readlines()
        header_lines = nps + 1
        for line in lines[header_lines:]:
            parts = line.strip().split()
            if len(parts) >= 3:  # Accept 3 or more columns (time, x, y, ...)
                try:
                    floats = [float(x) for x in parts]
                    data_lines.append(floats)
                except ValueError:
                    continue
    if not data_lines:
        print(f"Warning: No valid data found in {his_file}. Skipping.")
        return
    raw_data = np.array(data_lines)

    try:
        if raw_data.shape[1] == 3:
            data = raw_data.reshape(-1, nps, 3)
        elif raw_data.shape[1] == 4:
            data = raw_data.reshape(-1, nps, 4)
        elif raw_data.shape[1] == 5:
            data = raw_data.reshape(-1, nps, 5)
        else:
            print(f"Warning: Unexpected number of columns in {his_file}. Skipping.")
            return
    except ValueError as e:
        print(f"Error reshaping data for {his_file}: {e}. Check if the file is complete. Skipping.")
        return

    for prb in range(nps):
        print(f"Processing Probe {prb + 1}")
        t = data[:, prb, 0]  # Keep all data points
        # Extract both u and v components
        u = data[:, prb, 1] if data.shape[2] > 1 else data[:, prb, 0]
        v = data[:, prb, 2] if data.shape[2] > 2 else data[:, prb, 1]

        if len(t) < 2:
            print(f"Probe {prb + 1} has insufficient data. Skipping.")
            continue

        # Calculate means before interpolation
        u_mean, v_mean = np.mean(u), np.mean(v)
        u_zm, v_zm = u - u_mean, v - v_mean

        tn, un = interpolate_signal(t, u_zm)
        if len(tn) < 2:
            print(f"Probe {prb + 1} has insufficient data after interpolation. Skipping.")
            continue

        _, vn = interpolate_signal(t, v_zm, t_new=tn)
        dt = (tn[-1] - tn[0]) / (len(tn) - 1) if len(tn) > 1 else 0
        if dt == 0:
            print(f"Time step is zero for Probe {prb + 1}. Skipping.")
            continue

        fig = plt.figure(figsize=(3 * plot_params["fig_width"], plot_params["fig_height"]))
        gs = gridspec.GridSpec(2, 2, width_ratios=[2.2, 1])
        ax_signal = fig.add_subplot(gs[0, 0])
        ax_fft = fig.add_subplot(gs[1, 0])
        ax_phase = fig.add_subplot(gs[:, 1])

        ax_signal.scatter(tn, un, c="b", s=0.05, label=f"u (mean={u_mean:.5f})")
        ax_signal.scatter(tn, vn, c="r", s=0.05, label=f"v (mean={v_mean:.5f})")
        ax_signal.set_xlabel(r"$t$")
        ax_signal.set_ylabel(r"$u' / v'$")
        ax_signal.set_title(f"x,y,z={coords[prb][0]},{coords[prb][1]},{coords[prb][2]}")
        ax_signal.legend()

        plot_fft(ax_fft, tn, un, dt, label="u'", color="b", plot_reference=True)
        plot_fft(ax_fft, tn, vn, dt, label="v'", color="r")
        ax_fft.legend(loc="upper right")
        ax_fft.set_title("Spectrum")
        plot_phase_portrait(ax_phase, un, dt, color="b", name="u'", add_label=True)
        plot_phase_portrait(ax_phase, vn, dt, color="r", name="v'", add_label=True)
        ax_phase.legend()

        file_root = his_file.split(".")[0]
        fname = f"{file_root}_probe{prb + 1}_signals.{plot_params['format']}"
        print("Saving ", fname)
        fig.tight_layout()
        fig.savefig(fname, format=plot_params["format"], dpi=plot_params["quality"], bbox_inches=plot_params["adjust"])
        plt.close(fig)


def process_lift_file(lift_drag_file):
    """Load, process, and plot data from a lift_drag file."""
    print(f"--- Processing Lift/Drag File: {lift_drag_file} ---")
    loader = LiftDragLoader(lift_drag_file)
    t, dgx, dgy = loader.t, loader.dgx, loader.dgy

    if t.size < 2:
        print(f"No data loaded from {lift_drag_file}. Skipping.")
        return

    dgx_mean, dgy_mean = np.mean(dgx), np.mean(dgy)
    dgx_zm, dgy_zm = dgx - dgx_mean, dgy - dgy_mean

    tn, dgx_interp = interpolate_signal(t, dgx_zm)
    if len(tn) < 2:
        print(f"Insufficient data in {lift_drag_file} after interpolation. Skipping.")
        return

    dt = (tn[-1] - tn[0]) / (len(tn) - 1) if len(tn) > 1 else 0
    if dt == 0:
        print(f"Could not determine time step for {lift_drag_file}. Skipping.")
        return

    _, dgy_interp = interpolate_signal(t, dgy_zm, t_new=tn)

    fig = plt.figure(figsize=(3 * plot_params["fig_width"], plot_params["fig_height"]))
    gs = gridspec.GridSpec(2, 2, width_ratios=[2.2, 1])
    ax_signal = fig.add_subplot(gs[0, 0])
    ax_fft = fig.add_subplot(gs[1, 0])
    ax_phase = fig.add_subplot(gs[:, 1])

    ax_signal.scatter(tn, dgx_interp, c="b", s=0.05, label=f"Cx (mean={dgx_mean:.5f})")
    ax_signal.scatter(tn, dgy_interp, c="g", s=0.05, label=f"Cy (mean={dgy_mean:.5f})")
    ax_signal.set_xlabel(r"$t$")
    ax_signal.set_ylabel(r"$C_x' / C_y'$")
    ax_signal.set_title("Time vs Cx and Cy (zero-mean)")
    ax_signal.legend()

    plot_fft(ax_fft, tn, dgx_interp, dt, "Cx'", "b", scaling="spectrum", plot_reference=True)
    plot_fft(ax_fft, tn, dgy_interp, dt, "Cy'", "g", scaling="spectrum")
    ax_fft.legend(loc="upper right")
    ax_fft.set_title("Spectrum")

    plot_phase_portrait(ax_phase, dgx_interp, dt, "b", name="C_x'", add_label=True)
    plot_phase_portrait(ax_phase, dgy_interp, dt, "g", name="C_y'", add_label=True)
    ax_phase.legend()

    file_root = lift_drag_file.split(".")[0]
    fname = f"{file_root}_signals.{plot_params['format']}"
    print("Saving ", fname)
    fig.tight_layout()
    fig.savefig(fname, format=plot_params["format"], dpi=plot_params["quality"], bbox_inches=plot_params["adjust"])
    plt.close(fig)


# --- Main Execution ---


def main():
    parser = argparse.ArgumentParser(description="Process and plot signals from Nek5000 simulations.")
    parser.add_argument("--his", action="store_true", help="Process history point file (e.g., 1cyl.his).")
    parser.add_argument("--lift", action="store_true", help="Process lift/drag file (e.g., lift_drag.dat).")
    args = parser.parse_args()

    run_his = args.his
    run_lift = args.lift

    # If no flags are specified, try to run both
    if not run_his and not run_lift:
        print("No specific flag provided, attempting to process both his and lift files.")
        run_his = True
        run_lift = True

    if run_his:
        case_name = get_case_name_from_usr()
        his_file = f"{case_name}.all"
        if not os.path.isfile(his_file):
            his_file = f"{case_name}.his"
        if os.path.isfile(his_file):
            process_his_file(his_file)
        else:
            print("Info: No '.his' or '.all' file found, skipping history file processing.")

    if run_lift:
        lift_drag_file = "lift_drag.dat"
        if not os.path.isfile(lift_drag_file):
            lift_drag_file = "lift_drag.all"
        if os.path.isfile(lift_drag_file):
            process_lift_file(lift_drag_file)
        else:
            print("Info: No 'lift_drag.dat' or 'lift_drag.all' file found, skipping lift/drag processing.")


if __name__ == "__main__":
    main()
