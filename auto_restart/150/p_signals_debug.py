#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os, sys, argparse

try:
    sys.path.insert(0, os.path.join(os.environ["NEKSTAB_SOURCE_ROOT"], "bin"))
except KeyError:
    raise EnvironmentError("NEKSTAB_SOURCE_ROOT environment variable is not set.")
from nekStab_tools import (
    interpolate_signal,
    LiftDragLoader,
)

def plot_debug_signals_lift(lift_drag_file):
    """Debug plot for lift/drag showing time sampling analysis and interpolation quality."""
    print(f"--- Creating Debug Plot for: {lift_drag_file} ---")
    loader = LiftDragLoader(lift_drag_file)
    
    if len(loader.t) == 0:
        print("No lift_drag data found for debugging")
        return
        
    # Create debug plots with more useful info
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
    
    # Plot 1: Time step analysis
    if len(loader.t) > 1:
        dt_array = np.diff(loader.t)
        ax1.plot(loader.t[:-1], dt_array, 'ko-', markersize=2)
        ax1.axhline(y=np.mean(dt_array), color='r', linestyle='--', label=f'Mean dt = {np.mean(dt_array):.6f}')
        ax1.axhline(y=np.max(dt_array), color='orange', linestyle='--', label=f'Max dt = {np.max(dt_array):.6f}')
        ax1.axhline(y=np.min(dt_array), color='green', linestyle='--', label=f'Min dt = {np.min(dt_array):.6f}')
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Time Step (dt)')
        ax1.set_title('Time Step Analysis (FFT uses Max dt)')
        ax1.legend()
        ax1.grid(True)
    
    # Plot 2: Raw data quality check
    ax2.plot(loader.t, loader.dgx, 'b.-', label='Drag X', markersize=1, linewidth=0.5)
    ax2.plot(loader.t, loader.dgy, 'r.-', label='Lift Y', markersize=1, linewidth=0.5)
    # Mark potential discontinuities
    if len(loader.t) > 1:
        dt_array = np.diff(loader.t)
        large_gaps = np.where(dt_array > 2*np.mean(dt_array))[0]
        for gap_idx in large_gaps:
            ax2.axvline(x=loader.t[gap_idx], color='orange', linestyle=':', alpha=0.7, label='Potential gap' if gap_idx == large_gaps[0] else "")
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Force')
    ax2.set_title('Raw Data with Gap Detection')
    ax2.legend()
    ax2.grid(True)
    
    # Plot 3: FFT-safe vs naive interpolation comparison
    if len(loader.t) > 10:
        # Take middle 20% of data for detailed view
        start_idx = int(0.4 * len(loader.t))
        end_idx = int(0.6 * len(loader.t))
        t_zoom = loader.t[start_idx:end_idx]
        dgx_zoom = loader.dgx[start_idx:end_idx]
        
        # Compare FFT-safe (max dt) vs naive (high resolution) interpolation
        t_fft_safe, dgx_fft_safe = interpolate_signal(t_zoom, dgx_zoom, fft_safe=True)
        t_naive, dgx_naive = interpolate_signal(t_zoom, dgx_zoom, num_points=len(t_zoom)*5, fft_safe=False)
        
        ax3.plot(t_zoom, dgx_zoom, 'ko', label='Raw data', markersize=4)
        ax3.plot(t_fft_safe, dgx_fft_safe, 'b-', label='FFT-safe (max dt)', alpha=0.8, linewidth=2)
        ax3.plot(t_naive, dgx_naive, 'r--', label='Naive (fine grid)', alpha=0.6, linewidth=1)
        ax3.set_xlabel('Time')
        ax3.set_ylabel('Drag Force')
        ax3.set_title('FFT-Safe vs Naive Interpolation')
        ax3.legend()
        ax3.grid(True)
    
    # Plot 4: Signal statistics summary
    dt_array = np.diff(loader.t)
    dt_variability = np.std(dt_array) / np.mean(dt_array) * 100  # CV as percentage
    nyquist_max_dt = 0.5 / np.max(dt_array)  # Conservative Nyquist frequency
    
    stats_text = f"""Signal Statistics:
Length: {len(loader.t)} points
Time range: {loader.t[0]:.3f} to {loader.t[-1]:.3f}
Duration: {loader.t[-1] - loader.t[0]:.3f}

Forces:
Drag X: mean={np.mean(loader.dgx):.5f}, std={np.std(loader.dgx):.5f}
Lift Y: mean={np.mean(loader.dgy):.5f}, std={np.std(loader.dgy):.5f}

Time Step Analysis:
Mean dt: {np.mean(dt_array):.6f}
Min dt:  {np.min(dt_array):.6f}
Max dt:  {np.max(dt_array):.6f}
Std dt:  {np.std(dt_array):.6f}
Variability: {dt_variability:.1f}% (CV)

FFT Considerations:
Max resolvable freq (conservative): {nyquist_max_dt:.3f}
{'UNIFORM time step' if dt_variability < 1 else 'VARIABLE time step - use max dt for FFT'}"""
    
    ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes, fontsize=9,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.axis('off')
    ax4.set_title('Signal Statistics')
    
    plt.tight_layout()
    debug_fname = 'lift_drag_debug.png'
    plt.savefig(debug_fname, dpi=150, bbox_inches='tight')
    print(f"Debug plot saved as '{debug_fname}'")
    plt.close(fig)


def plot_debug_signals_his(his_file):
    """Debug plot for .his file showing probe data quality and interpolation."""
    print(f"--- Creating Debug Plot for: {his_file} ---")
    
    with open(his_file, "r") as file:
        try:
            nps = int(file.readline().rstrip())
        except (ValueError, IndexError):
            print(f"Error: Could not read number of probes from {his_file}.")
            return
        
        coords = np.zeros((nps, 3))
        for n, line in zip(range(nps), file):
            try:
                coord_parts = line.split()
                if len(coord_parts) >= 3:
                    coords[n] = coord_parts[:3]
                else:
                    coords[n] = [0, 0, 0]
            except (ValueError, IndexError):
                coords[n] = [0, 0, 0]

    data_lines = []
    with open(his_file, "r") as file:
        lines = file.readlines()
        header_lines = nps + 1
        for line in lines[header_lines:]:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    floats = [float(x) for x in parts]
                    data_lines.append(floats)
                except ValueError:
                    continue
    
    if not data_lines:
        print(f"Warning: No valid data found in {his_file}.")
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
            print(f"Warning: Unexpected number of columns in {his_file}.")
            return
    except ValueError as e:
        print(f"Error reshaping data for {his_file}: {e}.")
        return

    # Debug plot for first probe only
    prb = 0
    t = data[:, prb, 0]  # Keep all data points
    u = data[:, prb, 1] if data.shape[2] > 1 else data[:, prb, 0]
    v = data[:, prb, 2] if data.shape[2] > 2 else data[:, prb, 1]
    
    if len(t) < 2:
        print(f"Probe {prb + 1} has insufficient data.")
        return
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
    
    # Plot 1: Time step analysis
    dt_array = np.diff(t)
    ax1.plot(t[:-1], dt_array, 'ko-', markersize=2)
    ax1.axhline(y=np.mean(dt_array), color='r', linestyle='--', label=f'Mean dt = {np.mean(dt_array):.6f}')
    ax1.axhline(y=np.max(dt_array), color='orange', linestyle='--', label=f'Max dt = {np.max(dt_array):.6f}')
    ax1.axhline(y=np.min(dt_array), color='green', linestyle='--', label=f'Min dt = {np.min(dt_array):.6f}')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Time Step (dt)')
    ax1.set_title('Time Step Analysis (FFT uses Max dt)')
    ax1.legend()
    ax1.grid(True)
    
    # Plot 2: Raw velocity components
    ax2.plot(t, u, 'b.-', label='u velocity', markersize=1, linewidth=0.5)
    ax2.plot(t, v, 'r.-', label='v velocity', markersize=1, linewidth=0.5)
    # Mark potential discontinuities
    large_gaps = np.where(dt_array > 2*np.mean(dt_array))[0]
    for gap_idx in large_gaps:
        ax2.axvline(x=t[gap_idx], color='orange', linestyle=':', alpha=0.7, label='Potential gap' if gap_idx == large_gaps[0] else "")
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Velocity')
    ax2.set_title(f'Raw Velocity Data - Probe 1 at ({coords[0][0]:.3f}, {coords[0][1]:.3f}, {coords[0][2]:.3f})')
    ax2.legend()
    ax2.grid(True)
    
    # Plot 3: Zero-mean signals and interpolation check
    u_mean, v_mean = np.mean(u), np.mean(v)
    u_zm, v_zm = u - u_mean, v - v_mean
    
    # Take middle section for detailed view
    if len(t) > 20:
        start_idx = int(0.4 * len(t))
        end_idx = int(0.6 * len(t))
        t_zoom = t[start_idx:end_idx]
        u_zoom = u_zm[start_idx:end_idx]
        v_zoom = v_zm[start_idx:end_idx]
        
        t_interp, u_interp = interpolate_signal(t_zoom, u_zoom, num_points=len(t_zoom)*3)
        _, v_interp = interpolate_signal(t_zoom, v_zoom, num_points=len(t_zoom)*3)
        
        ax3.plot(t_zoom, u_zoom, 'bo', label='u\' (raw)', markersize=3)
        ax3.plot(t_interp, u_interp, 'b-', label='u\' (interp)', alpha=0.7)
        ax3.plot(t_zoom, v_zoom, 'ro', label='v\' (raw)', markersize=3)
        ax3.plot(t_interp, v_interp, 'r-', label='v\' (interp)', alpha=0.7)
        ax3.set_xlabel('Time')
        ax3.set_ylabel('Zero-mean Velocity')
        ax3.set_title('Zero-Mean Interpolation (Middle 20%)')
        ax3.legend()
        ax3.grid(True)
    
    # Plot 4: Statistics
    dt_variability = np.std(dt_array) / np.mean(dt_array) * 100  # CV as percentage
    nyquist_max_dt = 0.5 / np.max(dt_array)  # Conservative Nyquist frequency
    
    stats_text = f"""Probe 1 Statistics:
Position: ({coords[0][0]:.3f}, {coords[0][1]:.3f}, {coords[0][2]:.3f})
Length: {len(t)} points
Time range: {t[0]:.3f} to {t[-1]:.3f}
Duration: {t[-1] - t[0]:.3f}

Velocities:
u: mean={u_mean:.5f}, std={np.std(u):.5f}
v: mean={v_mean:.5f}, std={np.std(v):.5f}

Time Step Analysis:
Mean dt: {np.mean(dt_array):.6f}
Min dt:  {np.min(dt_array):.6f}
Max dt:  {np.max(dt_array):.6f}
Std dt:  {np.std(dt_array):.6f}
Variability: {dt_variability:.1f}% (CV)

FFT Considerations:
Max resolvable freq (conservative): {nyquist_max_dt:.3f}
{'UNIFORM time step' if dt_variability < 1 else 'VARIABLE time step - use max dt for FFT'}

Total probes: {nps}"""
    
    ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes, fontsize=9,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.axis('off')
    ax4.set_title('Probe Statistics')
    
    plt.tight_layout()
    debug_fname = f'{his_file.split(".")[0]}_debug.png'
    plt.savefig(debug_fname, dpi=150, bbox_inches='tight')
    print(f"Debug plot saved as '{debug_fname}'")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Debug plots for Nek5000 simulation signals.")
    parser.add_argument("--his", action="store_true", help="Create debug plot for history point file.")
    parser.add_argument("--lift", action="store_true", help="Create debug plot for lift/drag file.")
    args = parser.parse_args()

    if not args.his and not args.lift:
        print("Creating debug plots for both his and lift files.")
        args.his = True
        args.lift = True

    if args.lift:
        lift_drag_file = "lift_drag.all"
        if not os.path.isfile(lift_drag_file):
            lift_drag_file = "lift_drag.dat"
        if os.path.isfile(lift_drag_file):
            plot_debug_signals_lift(lift_drag_file)
        else:
            print("Info: No 'lift_drag.dat' or 'lift_drag.all' file found for debug plotting.")
    
    if args.his:
        his_file = "1cyl.all"
        if not os.path.isfile(his_file):
            his_file = "1cyl.his"
        if os.path.isfile(his_file):
            plot_debug_signals_his(his_file)
        else:
            print("Info: No '.his' or '.all' file found for debug plotting.")


if __name__ == "__main__":
    main()