#!/usr/bin/env python
from pickletools import string1
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
    read_his_file,
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
    'peak_detector_threshold': 0.01,  # Default threshold for peak detection in FFT
    'peak_matching_tolerance': 2.0,   # Tolerance for matching peaks to reference values (200%)
    # Spectrum plot configuration
    "xscale": "log",    # x-axis scale for spectrum plot ("log" or "linear")
    "yscale": "log",    # y-axis scale for spectrum plot ("log" or "linear")
    "show_peaks": True, # show peak markers and vertical lines
    "show_reference": True, # show reference Strouhal number lines
    "auto_xlim_from_peaks": True, # automatically set xlim_min based on smallest peak frequency
    "fft_scaling": "density", # FFT scaling: "spectrum" or "density"
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

# --- Reference Line Colors Configuration ---
# Define consistent colors for reference St values (can be easily changed here)
reference_colors = {
    "St_1": "red",      # St_1 is always red
    "St_2": "green",    # St_2 is always green  
    "St_3": "blue",     # St_3 is always blue
}

# --- Reference Values Dictionary ---
reference_values = {
    "50": {
        "St_1": 0.125, "description_1": "St_1 Re=50",
        "St_2": 0.2,   "description_2": "St_2 Re=50",
        "St_3": 0.3,   "description_3": "St_3 Re=50",
    },
    "150": {#340
        "St_1": 0.1316, "description_1": "St_1 Re=340",
        "St_2": 0.0377,    "description_2": "St_2 Re=340",
        "St_3": 0.0009,    "description_3": "St_3 Re=340",
    },
    "355": {
        "St_1": 0.1280, "description_1": "St_1 Re=355",
        "St_2": 0.0350, "description_2": "St_2 Re=355",
        "St_3": 0.0116, "description_3": "St_3 Re=355",
    },
    # Add more Reynolds numbers here as needed:
    # "100": {
    #     "St_1": 0.164, "description_1": "St_1 Re=100",
    #     "St_2": 0.21,  "description_2": "St_2 Re=100",
    # },
    # "150": {
    #     "St_1": 0.182, "description_1": "St_1 Re=150",
    # },
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


# --- Helper Functions ---
def find_closest_peaks_to_references(peak_freqs, peak_vals, ref_freqs, ref_labels, tolerance=0.1):
    """Find the closest data peaks to reference frequencies and return matching pairs."""
    matches = []
    
    for ref_freq, ref_label in zip(ref_freqs, ref_labels):
        if len(peak_freqs) == 0:
            continue
            
        # Find the closest peak to this reference frequency
        distances = [abs(peak - ref_freq) for peak in peak_freqs]
        min_distance_idx = np.argmin(distances)
        closest_peak_freq = peak_freqs[min_distance_idx]
        closest_peak_val = peak_vals[min_distance_idx]
        min_distance = distances[min_distance_idx]
        
        # Check if the closest peak is within tolerance (relative to reference frequency)
        relative_error = min_distance / ref_freq
        if relative_error <= tolerance:
            matches.append({
                'ref_freq': ref_freq,
                'ref_label': ref_label,
                'peak_freq': closest_peak_freq,
                'peak_val': closest_peak_val,
                'distance': min_distance,
                'relative_error': relative_error
            })
    
    return matches


# --- Plotting Functions ---
def plot_fft(ax, tn, vn, dt, label, color, plot_reference=False, case_reynolds="", set_limits=True):
    """Computes and plots the FFT (periodogram) of a signal."""
    scaling = signal_config.get("fft_scaling", "spectrum")
    ufreq, psd = periodogram_rfft(vn, fs=1.0 / dt, scaling=scaling)
    peak_freqs, peak_vals = find_peaks(ufreq, psd, threshold=signal_config.get('peak_detector_threshold', 0.01))

    # Plot reference line first so it appears first in legend (if available for this Re)
    show_reference = signal_config.get("show_reference", True)
    if plot_reference and show_reference and case_reynolds in reference_values:
        # Check if reference lines have already been plotted on this axis
        existing_labels = [line.get_label() for line in ax.get_lines()]
        ref_already_plotted = any(label.startswith("$St_") for label in existing_labels)
        
        if not ref_already_plotted:
            ref_data = reference_values[case_reynolds]
            # Only plot unique reference St values and avoid duplicate legend entries
            plotted_st = set()
            ref_freqs = []
            ref_labels = []
            
            for key in ref_data:
                if key.startswith("St_"):
                    st_num = key.split('_')[1]
                    st_value = ref_data[key]
                    if st_value in plotted_st:
                        continue
                    plotted_st.add(st_value)
                    st_key = f"St_{st_num}"
                    ref_color = reference_colors.get(st_key, "gray")  # Default to gray if not found
                    ref_label = f"$St_{{{st_num}}}={st_value}$"
                    ax.axvline(x=st_value, c=ref_color, lw=0.8, ls="--", label=ref_label)
                    ref_freqs.append(st_value)
                    ref_labels.append(st_key)
                    
            
            # Find closest peaks to reference values and mark them
            if len(peak_freqs) > 0 and len(ref_freqs) > 0:
                tolerance = signal_config.get('peak_matching_tolerance', 2.0)
                matches = find_closest_peaks_to_references(peak_freqs, peak_vals, ref_freqs, ref_labels, tolerance=tolerance)
                if matches:
                    print(f"Peak matching results: {len(matches)} matches out of {len(ref_freqs)} references")
                    for match in matches:
                        print(f"✓ {match['ref_label']}: Ref={match['ref_freq']:.4f}, Data={match['peak_freq']:.4f}, Δ={match['relative_error']*100:.1f}%")
                
                # Update reference labels to show data peak values for matches
                for match in matches:
                    # Find and update the corresponding reference line label
                    for line in ax.get_lines():
                        if line.get_label().startswith(f"$St_{{{match['ref_label'].split('_')[1]}}}"):
                            error_pct = match['relative_error'] * 100
                            new_label = f"$St_{{{match['ref_label'].split('_')[1]}}}={match['ref_freq']:.4f}$ ($\\Delta$={error_pct:.1f}%)"
                            line.set_label(new_label)
                            break

    lw = signal_config.get("line_width", 0.8)
    ls = signal_config.get("line_style", "-")
    mc = signal_config.get("fft_peak_marker_color", "k")
    xscale = signal_config.get("xscale", "log")
    yscale = signal_config.get("yscale", "log")
    
    # Plot based on scale configuration
    if xscale == "log" and yscale == "log":
        ax.loglog(ufreq, psd, c=color, lw=lw, ls=ls, label=label)
    elif xscale == "log":
        ax.semilogx(ufreq, psd, c=color, lw=lw, ls=ls, label=label)
    elif yscale == "log":
        ax.semilogy(ufreq, psd, c=color, lw=lw, ls=ls, label=label)
    else:
        ax.plot(ufreq, psd, c=color, lw=lw, ls=ls, label=label)
    
    show_peaks = signal_config.get("show_peaks", True)
    if show_peaks:
        ms = signal_config.get("fft_peak_marker_size", 40)
        ax.scatter(peak_freqs, peak_vals, c=mc, marker="x", s=ms)
        if len(peak_freqs) > 0:
            # Find the highest peak instead of the first peak
            highest_peak_idx = np.argmax(peak_vals)
            highest_peak_freq = peak_freqs[highest_peak_idx]
            ax.axvline(x=highest_peak_freq, c=color, lw=0.5, ls="--", label=f"$St={highest_peak_freq:.4f}$")

    # Only set limits if this is the last plot call or explicitly requested
    if set_limits:
        nyquist = 0.5 / dt
        xlim_min = signal_config.get("xlim_min") or 1e-2
        xlim_max = signal_config.get("xlim_max") or nyquist
        ax.set_xlim(xlim_min, xlim_max)
        
        ax.set_xlabel(r"$St$")
        if scaling == "spectrum":
            ax.set_ylabel("Power Spectrum")
        elif scaling == "density":
            ax.set_ylabel("PSD")
        else:
            ax.set_ylabel("Power")
    
    return psd, peak_freqs  # Return PSD and peak frequencies for limit calculation


def set_fft_limits(ax, all_psds, dt, all_peak_freqs=None, case_reynolds=None):
    """Set x and y-axis limits for FFT plot based on all PSDs, peak frequencies, and reference peaks."""
    nyquist = 0.5 / dt
    xlim_min = signal_config.get("xlim_min")
    xlim_max = signal_config.get("xlim_max") or nyquist
    auto_xlim = signal_config.get("auto_xlim_from_peaks", True)

    if xlim_min is None and auto_xlim:
        peak_freqs = [freq for peaks in all_peak_freqs for freq in peaks if freq > 0] if all_peak_freqs else []
        ref_freqs = []
        if case_reynolds is not None and case_reynolds in reference_values:
            ref_data = reference_values[case_reynolds]
            for key in ref_data:
                if key.startswith("St_"):
                    ref_freqs.append(ref_data[key])
        
        # Find the minimum between smallest peak and smallest reference marker
        min_peak = min(peak_freqs) if peak_freqs else float('inf')
        min_ref = min(ref_freqs) if ref_freqs else float('inf')
        
        if min_peak != float('inf') or min_ref != float('inf'):
            min_freq = min(min_peak, min_ref)
            xlim_min = min_freq / 10  # One decade to the left of the smallest frequency
        else:
            xlim_min = 1e-2
    else:
        xlim_min = xlim_min or 1e-2

    ax.set_xlim(xlim_min, xlim_max)
    
    # Set y-axis limits
    ylim_min = signal_config.get("ylim_min")
    ylim_max = signal_config.get("ylim_max")
    
    if not (ylim_min and ylim_max) and all_psds:
        # Combine all PSDs and find global min/max
        all_psd_vals = []
        for psd in all_psds:
            if len(psd) > 4:
                all_psd_vals.extend(psd[4:])  # Skip first few points
        
        if all_psd_vals:
            global_min = min(all_psd_vals)
            global_max = max(all_psd_vals) * 10
            auto_min = ylim_min or global_min
            auto_max = ylim_max or global_max
            ax.set_ylim(auto_min, auto_max)


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
    
    # Read the file using the new function
    coords, t, data = read_his_file(his_file)
    if coords is None:
        return
    
    nps = coords.shape[0]
    n_fields = data.shape[2]
    if3d = n_fields > 2

    # Determine which probes to process
    probes_to_process = range(nps) if plot_all_probes else [0]
    print(f"Processing {'all' if plot_all_probes else 'first'} probe{'s' if plot_all_probes else ''}")
    
    # Get config values once
    tmin, tmax = signal_config.get("tmin"), signal_config.get("tmax")
    marker_size = signal_config.get("marker_size", 0.05)
    marker_color_u = signal_config.get("marker_color_u", "b")
    marker_color_v = signal_config.get("marker_color_v", "r")

    for prb in probes_to_process:
        # Extract velocity components
        u = data[:, prb, 0] if n_fields > 0 else np.zeros_like(t)
        v = data[:, prb, 1] if n_fields > 1 else np.zeros_like(t)
        
        if if3d and n_fields > 2 and flip:
            w = data[:, prb, 2]
            v, w = w, v

        # Skip if insufficient data
        if len(t) < 2:
            print(f"Probe {prb + 1} has insufficient data. Skipping.")
            continue

        # Calculate means and remove them
        u_mean, v_mean = np.mean(u), np.mean(v)
        u_zm, v_zm = u - u_mean, v - v_mean

        # Interpolate signals
        tn, un = interpolate_signal(t, u_zm, tmin=tmin, tmax=tmax)

        _, vn = interpolate_signal(t, v_zm, t_new=tn)
        dt = (tn[-1] - tn[0]) / (len(tn) - 1)
        if dt == 0:
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

        # Plot FFTs and collect PSDs and peak frequencies for proper limit setting
        psd_u, peaks_u = plot_fft(ax_fft, tn, un, dt, label="u'", color=marker_color_u, plot_reference=True, case_reynolds=case_reynolds, set_limits=False)
        psd_v, peaks_v = plot_fft(ax_fft, tn, vn, dt, label="v'", color=marker_color_v, case_reynolds=case_reynolds, set_limits=False)
        
        # Set limits based on all PSDs and peak frequencies
        set_fft_limits(ax_fft, [psd_u, psd_v], dt, [peaks_u, peaks_v], case_reynolds)
        ax_fft.set_xlabel(r"$St$")
        scaling = signal_config.get("fft_scaling", "spectrum")
        if scaling == "spectrum":
            ax_fft.set_ylabel("Power Spectrum")
        elif scaling == "density":
            ax_fft.set_ylabel("PSD")
        else:
            ax_fft.set_ylabel("Power")
        ax_fft.legend(loc="lower right")
        ax_fft.set_title("Spectrum")
        plot_phase_portrait(ax_phase, un, dt, color=marker_color_u, name="u'", add_label=True)
        plot_phase_portrait(ax_phase, vn, dt, color=marker_color_v, name="v'", add_label=True)
        ax_phase.set_aspect('equal')  # Set equal aspect ratio after both plots
        ax_phase.legend()
        # Let matplotlib autoscale phase portrait axes

        # Save figure
        fname = f"{his_file.split('.')[0]}_probe{prb + 1}_signals.{plot_params['format']}"
        fig.tight_layout()
        fig.savefig(fname, format=plot_params["format"], dpi=plot_params["quality"], bbox_inches=plot_params["adjust"])
        plt.close(fig)


def process_lift_file(lift_drag_file, case_reynolds="", flip=False):
    """Load, process, and plot data from a lift_drag file."""
    print(f"--- Processing: {lift_drag_file} ---")
    loader = LiftDragLoader(lift_drag_file, flip=flip)
    t, dgx, dgy = loader.t, loader.dgx, loader.dgy

    if t.size < 2:
        print(f"No data loaded from {lift_drag_file}. Skipping.")
        return

    # Calculate means and remove them
    dgx_mean, dgy_mean = np.mean(dgx), np.mean(dgy)
    dgx_zm, dgy_zm = dgx - dgx_mean, dgy - dgy_mean

    # Interpolate signals
    tmin, tmax = signal_config.get("tmin"), signal_config.get("tmax")
    tn, dgx_interp = interpolate_signal(t, dgx_zm, tmin=tmin, tmax=tmax)
    if len(tn) < 2:
        return

    dt = (tn[-1] - tn[0]) / (len(tn) - 1)
    if dt == 0:
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

    # Plot FFTs and collect PSDs and peak frequencies for proper limit setting
    psd_cx, peaks_cx = plot_fft(ax_fft, tn, dgx_interp, dt, "Cx'", "b", plot_reference=True, case_reynolds=case_reynolds, set_limits=False)
    psd_cy, peaks_cy = plot_fft(ax_fft, tn, dgy_interp, dt, "Cy'", "r", plot_reference=False, case_reynolds=case_reynolds, set_limits=False)
    
    # Set limits based on all PSDs and peak frequencies
    set_fft_limits(ax_fft, [psd_cx, psd_cy], dt, [peaks_cx, peaks_cy], case_reynolds)
    ax_fft.set_xlabel(r"$St$")
    scaling = signal_config.get("fft_scaling", "spectrum")
    if scaling == "spectrum":
        ax_fft.set_ylabel("Power Spectrum")
    elif scaling == "density":
        ax_fft.set_ylabel("PSD")
    else:
        ax_fft.set_ylabel("Power")
    ax_fft.legend(loc="lower right")
    ax_fft.set_title("Spectrum")

    plot_phase_portrait(ax_phase, dgx_interp, dt, "b", name="C_x'", add_label=True)
    plot_phase_portrait(ax_phase, dgy_interp, dt, "r", name="C_y'", add_label=True)
    ax_phase.set_aspect('equal')  # Set equal aspect ratio after both plots
    ax_phase.legend()

    # Save figure
    fname = f"{lift_drag_file.split('.')[0]}_signals.{plot_params['format']}"
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
