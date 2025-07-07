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

# --- Signal/Plot Configuration ---
signal_config = {
    "tmin": None,         # minimum time for cropping (None = auto)
    "tmax": None,       # maximum time for cropping (None = auto)
    "xlim_min": None,   # min x-axis limit for plots (None = auto)
    "xlim_max": None,   # max x-axis limit for plots (None = auto)
    "ylim_min": None,   # min y-axis limit for plots (None = auto)
    "ylim_max": None,   # max y-axis limit for plots (None = auto)
    "marker_size": 0.05, # marker size for scatter plots
    "marker_color_u": "b", # color for u/Cx
    "marker_color_v": "r", # color for v/Cy
    "line_width": 0.8,  # line width for lines
    "line_style": "-", # line style for lines
    "fft_peak_marker_size": 40, # marker size for FFT peaks
    'fft_peak_marker_color': 'red',
    'peak_detector_threshold': 0.01,  # Default threshold for peak detection in FFT # color for FFT peaks
}

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

# --- Reference Values Dictionary ---
reference_values = {
    "50": {"strouhal": 0.125, "description": "Re=50"},
    # Add more Reynolds numbers here as needed:
    # "100": {"strouhal": 0.164, "description": "Re=100"},
    # "150": {"strouhal": 0.182, "description": "Re=150"},
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
def plot_fft(ax, tn, vn, dt, label, color, scaling="spectrum", plot_reference=False, case_reynolds=""):
    """Computes and plots the FFT (periodogram) of a signal."""
    ufreq, psd = periodogram_rfft(vn, fs=1.0 / dt, scaling=scaling)
    peak_freqs, peak_vals = find_peaks(ufreq, psd, threshold=signal_config.get('peak_detector_threshold', 0.01))

    # Plot reference line first so it appears first in legend (if available for this Re)
    if plot_reference and case_reynolds in reference_values:
        ref_data = reference_values[case_reynolds]
        strouhal = ref_data["strouhal"]
        description = ref_data["description"]
        ax.axvline(x=strouhal, c="k", lw=0.5, ls="--", label=f"Ref. ${description}$ $St={strouhal}$")

    lw = signal_config.get("line_width", 0.8)
    ls = signal_config.get("line_style", "-")
    ms = signal_config.get("fft_peak_marker_size", 40)
    mc = signal_config.get("fft_peak_marker_color", "k")
    ax.loglog(ufreq, psd, c=color, lw=lw, ls=ls, label=label)
    ax.scatter(peak_freqs, peak_vals, c=mc, marker="x", s=ms)
    if len(peak_freqs) > 0:
        # Find the highest peak instead of the first peak
        highest_peak_idx = np.argmax(peak_vals)
        highest_peak_freq = peak_freqs[highest_peak_idx]
        ax.axvline(x=highest_peak_freq, c=color, lw=0.5, ls="--", label=f"$St={highest_peak_freq:.4f}$")

    nyquist = 0.5 / dt
    xlim_min = signal_config.get("xlim_min", 1e-2)
    xlim_max = signal_config.get("xlim_max", nyquist)
    ax.set_xlim(xlim_min if xlim_min is not None else 1e-2, xlim_max if xlim_max is not None else nyquist)
    psd_without_0 = psd[4:]
    ylim_min = signal_config.get("ylim_min", None)
    ylim_max = signal_config.get("ylim_max", None)
    if len(psd_without_0) > 0:
        auto_min = psd_without_0.min()
        auto_max = psd_without_0.max() * 10
        ax.set_ylim(ylim_min if ylim_min is not None else auto_min, ylim_max if ylim_max is not None else auto_max)
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

    # Add reference lines at origin
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    lw = signal_config.get("line_width", 0.8)
    ls = signal_config.get("line_style", "-")
    ms = signal_config.get("marker_size", 0.05)
    ax.set_xlabel(f"${name}$")
    ax.set_ylabel(rf"$\dot{{{name}}}$")
    ax.set_title("Phase Portrait")
    if add_label:
        ax.plot(v_trim, vdot_trim, c=color, lw=lw, ls=ls, label=f"${name}$")
    else:
        ax.plot(v_trim, vdot_trim, c=color, lw=lw, ls=ls)


# --- Data Processors ---
def process_his_file(his_file, plot_all_probes=False, case_reynolds="", flip=False):
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
            if3d = False
        elif raw_data.shape[1] == 4:
            data = raw_data.reshape(-1, nps, 4)
            if3d = True
        elif raw_data.shape[1] == 5:
            data = raw_data.reshape(-1, nps, 5)
            if3d = True
        else:
            print(f"Warning: Unexpected number of columns in {his_file}. Skipping.")
            return
    except ValueError as e:
        print(f"Error reshaping data for {his_file}: {e}. Check if the file is complete. Skipping.")
        return

    # Determine which probes to process
    if plot_all_probes:
        probes_to_process = range(nps)
        print(f"Processing all {nps} probes")
    else:
        probes_to_process = [0]  # Only first probe
        print(f"Processing only first probe (use --all to process all {nps} probes)")
    
    tmin = signal_config.get("tmin", None)
    tmax = signal_config.get("tmax", None)
    marker_size = signal_config.get("marker_size", 0.05)
    marker_color_u = signal_config.get("marker_color_u", "b")
    marker_color_v = signal_config.get("marker_color_v", "r")
    line_width = signal_config.get("line_width", 0.8)
    line_style = signal_config.get("line_style", "-")

    for prb in probes_to_process:
        print(f"Processing Probe {prb + 1}")
        t = data[:, prb, 0]  # Keep all data points
        # Extract both u and v components
        u = data[:, prb, 1] if data.shape[2] > 1 else data[:, prb, 0]
        v = data[:, prb, 2] if data.shape[2] > 2 else data[:, prb, 1]
        if if3d:
            w = data[:, prb, 3] if data.shape[2] > 3 else data[:, prb, 2]
            if flip:
                v, w = w, v 

        # Do NOT crop signals here! Cropping is handled in interpolate_signal.
        if len(t) < 2:
            print(f"Probe {prb + 1} has insufficient data. Skipping.")
            continue

        # Calculate means before interpolation
        u_mean, v_mean = np.mean(u), np.mean(v)
        u_zm, v_zm = u - u_mean, v - v_mean

        tn, un = interpolate_signal(t, u_zm, tmin=tmin, tmax=tmax)
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

        ax_signal.scatter(tn, un, c=marker_color_u, s=marker_size, label=f"u (mean={u_mean:.5f})")
        ax_signal.scatter(tn, vn, c=marker_color_v, s=marker_size, label=f"v (mean={v_mean:.5f})")
        ax_signal.set_xlabel(r"$t$")
        ax_signal.set_ylabel(r"$u' / v'$")
        ax_signal.set_title(f"Re={case_reynolds}, probe at x,y,z={coords[prb][0]},{coords[prb][1]},{coords[prb][2]}")
        ax_signal.legend()
        # For the time series plot, set axis limits automatically based on the signal
        ax_signal.set_xlim(tn[0], tn[-1])
        # Let matplotlib autoscale y-axis for the signal plot

        plot_fft(ax_fft, tn, un, dt, label="u'", color=marker_color_u, plot_reference=True, case_reynolds=case_reynolds)
        plot_fft(ax_fft, tn, vn, dt, label="v'", color=marker_color_v, case_reynolds=case_reynolds)
        ax_fft.legend(loc="upper right")
        ax_fft.set_title("Spectrum")
        plot_phase_portrait(ax_phase, un, dt, color=marker_color_u, name="u'", add_label=True)
        plot_phase_portrait(ax_phase, vn, dt, color=marker_color_v, name="v'", add_label=True)
        ax_phase.set_aspect('equal')  # Set equal aspect ratio after both plots
        ax_phase.legend()
        # xlim/ylim for phase portrait
        if xlim_min is not None or xlim_max is not None:
            ax_phase.set_xlim(left=xlim_min if xlim_min is not None else None, right=xlim_max if xlim_max is not None else None)
        if ylim_min is not None or ylim_max is not None:
            ax_phase.set_ylim(bottom=ylim_min if ylim_min is not None else None, top=ylim_max if ylim_max is not None else None)

        file_root = his_file.split(".")[0]
        fname = f"{file_root}_probe{prb + 1}_signals.{plot_params['format']}"
        print("Saving ", fname)
        fig.tight_layout()
        fig.savefig(fname, format=plot_params["format"], dpi=plot_params["quality"], bbox_inches=plot_params["adjust"])
        plt.close(fig)


def process_lift_file(lift_drag_file, case_reynolds="", flip=False):
    """Load, process, and plot data from a lift_drag file."""
    print(f"--- Processing Lift/Drag File: {lift_drag_file} ---")
    loader = LiftDragLoader(lift_drag_file, flip=flip)
    t, dgx, dgy = loader.t, loader.dgx, loader.dgy

    if t.size < 2:
        print(f"No data loaded from {lift_drag_file}. Skipping.")
        return

    dgx_mean, dgy_mean = np.mean(dgx), np.mean(dgy)
    dgx_zm, dgy_zm = dgx - dgx_mean, dgy - dgy_mean

    tmin = signal_config.get("tmin", None)
    tmax = signal_config.get("tmax", None)
    tn, dgx_interp = interpolate_signal(t, dgx_zm, tmin=tmin, tmax=tmax)
    if len(tn) < 2:
        print(f"Insufficient data in {lift_drag_file} after interpolation. Skipping.")
        return

    dt = (tn[-1] - tn[0]) / (len(tn) - 1) if len(tn) > 1 else 0
    if dt == 0:
        print(f"Could not determine time step for {lift_drag_file}. Skipping.")
        return

    _, dgy_interp = interpolate_signal(t, dgy_zm, t_new=tn, tmin=tmin, tmax=tmax)

    fig = plt.figure(figsize=(3 * plot_params["fig_width"], plot_params["fig_height"]))
    gs = gridspec.GridSpec(2, 2, width_ratios=[2.2, 1])
    ax_signal = fig.add_subplot(gs[0, 0])
    ax_fft = fig.add_subplot(gs[1, 0])
    ax_phase = fig.add_subplot(gs[:, 1])

    ax_signal.scatter(tn, dgx_interp, c="b", s=0.05, label=f"Cx (mean={dgx_mean:.5f})")
    ax_signal.scatter(tn, dgy_interp, c="r", s=0.05, label=f"Cy (mean={dgy_mean:.5f})")
    ax_signal.set_xlabel(r"$t$")
    ax_signal.set_ylabel(r"$C_x' / C_y'$")
    ax_signal.set_title(f"Re={case_reynolds}, Time vs Cx and Cy (zero-mean)")
    ax_signal.legend()

    plot_fft(ax_fft, tn, dgx_interp, dt, "Cx'", "b", scaling="spectrum", plot_reference=True, case_reynolds=case_reynolds)
    plot_fft(ax_fft, tn, dgy_interp, dt, "Cy'", "r", scaling="spectrum", case_reynolds=case_reynolds)
    ax_fft.legend(loc="upper right")
    ax_fft.set_title("Spectrum")

    plot_phase_portrait(ax_phase, dgx_interp, dt, "b", name="C_x'", add_label=True)
    plot_phase_portrait(ax_phase, dgy_interp, dt, "r", name="C_y'", add_label=True)
    ax_phase.set_aspect('equal')  # Set equal aspect ratio after both plots
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
    parser.add_argument("--all", action="store_true", help="Process all probes (default: only first probe).")
    args = parser.parse_args()

    run_his = args.his
    run_lift = args.lift
    plot_all_probes = args.all
    
    # Get Reynolds number from directory name (same as check_restart.py)
    working_directory = os.getcwd()
    case_reynolds = os.path.basename(working_directory)

    # If no flags are specified, try to run both
    if not run_his and not run_lift:
        print("No specific flag provided, attempting to process both his and lift files.")
        run_his = True
        run_lift = True

    case_name = get_case_name_from_usr()
    flip = (case_name.lower() == "sphere")

    if run_his:
        his_file = f"{case_name}.all"
        if not os.path.isfile(his_file):
            his_file = f"{case_name}.his"
        if os.path.isfile(his_file):
            process_his_file(his_file, plot_all_probes, case_reynolds, flip=flip)
        else:
            print("Info: No '.his' or '.all' file found, skipping history file processing.")

    if run_lift:
        lift_drag_file = "lift_drag.all"
        if not os.path.isfile(lift_drag_file):
            lift_drag_file = "lift_drag.dat"
        if os.path.isfile(lift_drag_file):
            process_lift_file(lift_drag_file, case_reynolds, flip=flip)
        else:
            print("Info: No 'lift_drag.dat' or 'lift_drag.all' file found, skipping lift/drag processing.")


if __name__ == "__main__":
    main()
