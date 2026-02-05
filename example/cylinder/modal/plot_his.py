#!/usr/bin/env python3
"""
p_his.py: Parse and plot Nek5000 history point (.his) files.

OUTPUTS: <basename><probe_num>_fft.png, <basename><probe_num>_phase_space.png
USAGE:   python p_his.py                     # Default: 1cyl.his
         python p_his.py file.his            # Specific file
         python p_his.py file.his --tskip 50 # Skip t < 50
         python p_his.py file.his --tmax 200 # Crop t > 200
         python p_his.py *.his               # Process multiple files
"""
import argparse
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d


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


def parse_his_file(filepath):
    """
    Parse a Nek5000 .his file.

    Returns:
        coords: (nps, ndim) array of probe coordinates
        data: (nsamples, nps, ncols) array where ncols is 4 (2D) or 5 (3D)
    """
    filepath = Path(filepath)

    with open(filepath, "r") as f:
        # First line: number of probes
        nps = int(f.readline().strip())

        # Next nps lines: probe coordinates (2 or 3 values)
        coords = []
        for _ in range(nps):
            vals = list(map(float, f.readline().split()))
            coords.append(vals)
        coords = np.array(coords)

    # Read data section with retry logic (for race condition with running sim)
    max_retries = 3
    for attempt in range(max_retries):
        try:
            raw_data = np.loadtxt(filepath, skiprows=nps + 1)

            # Determine data columns: 4 for 2D (t,u,v,p), 5 for 3D (t,u,v,w,p)
            if nps == 1:
                ncols = raw_data.shape[1]
                data = raw_data.reshape(-1, 1, ncols)
            else:
                ncols = raw_data.shape[1] // nps if raw_data.ndim > 1 else 4
                data = raw_data.reshape(-1, nps, ncols)

            # Validate: time should be non-negative
            t = data[1:, 0, 0]  # Skip first sample
            if t.size > 0 and t.min() >= 0:
                return coords, data

        except (ValueError, IndexError):
            pass

        if attempt < max_retries - 1:
            time.sleep(0.1)

    raise RuntimeError(f"Could not read valid data from {filepath} after {max_retries} attempts")


def format_coords(coord):
    """Format probe coordinates as string for display."""
    if len(coord) == 2:
        return f"x={coord[0]:.2f}, y={coord[1]:.2f}"
    return f"x={coord[0]:.2f}, y={coord[1]:.2f}, z={coord[2]:.2f}"


def process_probe(data, prb, coords, tskip=0, tmax=0):
    """
    Process a single probe's data and generate plots.

    Returns:
        dict with time series info and dominant frequency
    """
    # Extract time and v-velocity (skip first row)
    t = data[1:, prb, 0]
    v = data[1:, prb, 2]

    # Time range filtering
    tmin_orig, tmax_orig = t.min(), t.max()
    tmin = max(tmin_orig, tskip)
    tmax_use = min(tmax_orig, tmax) if tmax > 0 else tmax_orig

    print(f"  Time range: [{tmin_orig:.2f}, {tmax_orig:.2f}] -> [{tmin:.2f}, {tmax_use:.2f}]")

    # Crop indices
    itmin = np.searchsorted(t, tmin)
    itmax = np.searchsorted(t, tmax_use)
    if itmax <= itmin:
        print(f"  Warning: No data in time range, skipping")
        return None

    to = t[itmin:itmax]
    vo = v[itmin:itmax]

    # Compute uniform time step (use max dt for safety)
    dt_arr = np.diff(to)
    dt = dt_arr.max() if len(dt_arr) > 0 else 1.0
    tn = np.linspace(to[0], to[-1], int((to[-1] - to[0]) / dt), endpoint=False)

    # Remove mean and interpolate to uniform grid
    vo_mean = np.mean(vo)
    vo_centered = vo - vo_mean
    vn = interp1d(to, vo_centered, kind="nearest", fill_value="extrapolate")(tn)

    print(f"  v: mean={vo_mean:.6f}, range=[{vo_centered.min():.6f}, {vo_centered.max():.6f}]")

    # Compute PSD
    freqs, psd = periodogram_rfft(vn, fs=1.0/dt)
    peak_freqs, peak_vals = find_peaks(freqs, psd, threshold=0.01)

    dominant_st = peak_freqs[0] if len(peak_freqs) > 0 else None
    if dominant_st:
        print(f"  Dominant St = {dominant_st:.4f}")
        print(f"  St = {dominant_st}")  # Full precision for copy-paste

    return {
        "tn": tn,
        "vn": vn,
        "freqs": freqs,
        "psd": psd,
        "peak_freqs": peak_freqs,
        "peak_vals": peak_vals,
        "dominant_st": dominant_st,
        "dt": dt,
        "v_full": v,
        "coords": coords[prb],
    }


def plot_fft(result, output_path, ref_st=None):
    """Generate time series and PSD plot."""
    fig, axs = plt.subplots(2, figsize=(10, 6))

    # Time series
    axs[0].scatter(result["tn"], result["vn"], c="r", s=0.5)
    axs[0].set_xlabel(r"$t$")
    axs[0].set_ylabel(r"$v - \bar{v}$")
    axs[0].set_title(format_coords(result["coords"]))

    # PSD
    freqs, psd = result["freqs"], result["psd"]
    axs[1].semilogy(freqs, psd, c="r", lw=0.6)
    axs[1].scatter(result["peak_freqs"], result["peak_vals"], c="k", marker="x", s=40, label="peaks")

    if result["dominant_st"]:
        axs[1].axvline(x=result["dominant_st"], c="r", lw=0.5, ls="--",
                       label=f"$St={result['dominant_st']:.4f}$")

    if ref_st:
        axs[1].axvline(x=ref_st, c="k", lw=0.5, ls="--", label=f"Ref. $St={ref_st}$")

    axs[1].set_xlabel(r"$St$")
    axs[1].set_xlim(1e-4, 1)

    # Set y-limits, avoiding the zero-frequency spike
    psd_trim = psd[4:] if len(psd) > 4 else psd
    if len(psd_trim) > 0 and psd_trim.max() > 0:
        axs[1].set_ylim(psd_trim.min() * 0.1, psd_trim.max() * 10)

    axs[1].legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_phase_space(result, output_path):
    """Generate phase space plot (v vs dv/dt)."""
    v = result["v_full"]
    dt = result["dt"]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(v[2:-2], np.gradient(v)[2:-2] / dt, c="r", lw=0.8)
    ax.set_xlabel(r"$v$")
    ax.set_ylabel(r"$\dot{v}$")
    ax.set_title(f"Phase space: {format_coords(result['coords'])}")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Parse and plot Nek5000 history point (.his) files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python p_his.py                        # Process 1cyl.his
  python p_his.py myfile.his             # Process specific file
  python p_his.py *.his                  # Process all .his files
  python p_his.py file.his --tskip 100   # Skip t < 100
  python p_his.py file.his --tmax 500    # Crop t > 500
  python p_his.py file.his --ref-st 0.2  # Add reference Strouhal line
        """
    )
    parser.add_argument("files", nargs="*", default=["1cyl.his"],
                        help="Input .his file(s) (default: 1cyl.his)")
    parser.add_argument("--tskip", type=float, default=0,
                        help="Skip time values below this threshold")
    parser.add_argument("--tmax", type=float, default=0,
                        help="Crop time values above this threshold (0 = no limit)")
    parser.add_argument("--ref-st", type=float, default=None,
                        help="Reference Strouhal number to mark on PSD plot")
    parser.add_argument("--probes", type=str, default=None,
                        help="Probe indices to process, e.g., '1,3,5' or '1-5' (1-indexed, default: all)")

    args = parser.parse_args()

    for filepath in args.files:
        filepath = Path(filepath)
        if not filepath.exists():
            print(f"File not found: {filepath}")
            continue

        print(f"\nProcessing: {filepath}")

        try:
            coords, data = parse_his_file(filepath)
        except RuntimeError as e:
            print(f"  Error: {e}")
            continue

        nps = coords.shape[0]
        nsamples = data.shape[0]
        print(f"  Probes: {nps}, Dimensions: {coords.shape[1]}D, Samples: {nsamples}")

        # Parse probe selection
        if args.probes:
            probe_indices = []
            for part in args.probes.split(","):
                if "-" in part:
                    start, end = map(int, part.split("-"))
                    probe_indices.extend(range(start - 1, end))  # Convert to 0-indexed
                else:
                    probe_indices.append(int(part) - 1)
            probe_indices = [i for i in probe_indices if 0 <= i < nps]
        else:
            probe_indices = range(nps)

        basename = filepath.stem

        for prb in probe_indices:
            print(f"\n  Probe {prb + 1}/{nps}:")

            result = process_probe(data, prb, coords,
                                   tskip=args.tskip, tmax=args.tmax)
            if result is None:
                continue

            # Generate plots
            fft_path = filepath.parent / f"{basename}{prb + 1}_fft.png"
            phase_path = filepath.parent / f"{basename}{prb + 1}_phase_space.png"

            plot_fft(result, fft_path, ref_st=args.ref_st)
            plot_phase_space(result, phase_path)


if __name__ == "__main__":
    main()
