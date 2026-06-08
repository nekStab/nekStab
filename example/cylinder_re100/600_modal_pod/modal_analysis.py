#!/usr/bin/env python
"""
Modal analysis script: POD, DMD, SPOD on cylinder Re=100 wake snapshots.
Generates gallery-ready PNG plots (600 dpi, square-ish field panels).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
import tempfile

from pymodal import PODAnalyzer, DMDAnalyzer, SPODAnalyzer

# Configuration
CASE_DIR = Path(__file__).parent
NPZ_PATH = CASE_DIR / 'modal_snapshots.npz'
N_MODES = 10
SPOD_NFFT = 64
SPOD_OVERLAP = 0.5
DPI = 600
CYLINDER_RADIUS = 0.5

# Plot parameters
XLIM = (-2, 20)
YLIM = (-4, 4)


def load_grid():
    """Load grid from NPZ."""
    data = np.load(NPZ_PATH)
    x = data['x']
    y = data['y']
    data.close()
    return x, y


def reshape_mode(mode, nx, ny):
    """Reshape flattened mode (nx*ny,) to spatial grid (nx, ny)."""
    return mode.reshape(nx, ny)


def plot_field(ax, x, y, field, title, cmap = 'RdBu_r'):
    """Plot a 2D field with cylinder patch."""
    # field is (nx, ny); transpose for pcolormesh so x is horizontal
    im = ax.pcolormesh(x, y, field.T, cmap = cmap, shading = 'auto')
    
    # Add cylinder patch at origin (radius = 0.5)
    circle = patches.Circle((0, 0), CYLINDER_RADIUS, fill = True, 
                             color = 'white', edgecolor = 'black', linewidth = 0.5)
    ax.add_patch(circle)
    
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(title, fontsize = 10)
    plt.colorbar(im, ax = ax, label = 'magnitude')
    return im


def run_pod(npz_path, temp_dir):
    """Run POD and return analyzer."""
    print("Running POD analysis...")
    analyzer = PODAnalyzer(
        npz_path,
        results_dir = str(temp_dir),
        figures_dir = str(temp_dir),
        n_modes_save = N_MODES,
        use_parallel = False
    )
    analyzer.load_and_preprocess()
    analyzer.perform_pod()
    print(f"  Computed {analyzer.eigenvalues.shape[0]} modes")
    return analyzer


def run_dmd(npz_path, temp_dir):
    """Run DMD and return analyzer."""
    print("Running DMD analysis...")
    analyzer = DMDAnalyzer(
        npz_path,
        results_dir = str(temp_dir),
        figures_dir = str(temp_dir),
        n_modes_save = N_MODES,
        use_parallel = False
    )
    analyzer.load_and_preprocess()
    analyzer.perform_dmd()
    print(f"  Computed {analyzer.eigenvalues.shape[0]} modes")
    return analyzer


def run_spod(npz_path, temp_dir):
    """Run SPOD and return analyzer."""
    print("Running SPOD analysis...")
    analyzer = SPODAnalyzer(
        npz_path,
        nfft = SPOD_NFFT,
        overlap = SPOD_OVERLAP,
        results_dir = str(temp_dir),
        figures_dir = str(temp_dir),
        use_parallel = False
    )
    analyzer.load_and_preprocess()
    analyzer.run()  # Compute FFT blocks
    analyzer.perform_spod()
    print(f"  Computed modes shape {analyzer.modes.shape}")
    return analyzer


def plot_pod_spectrum(analyzer, x, y, output_dir):
    """Plot POD energy spectrum."""
    fig, ax = plt.subplots(figsize = (7, 6))
    
    eigenvalues = analyzer.eigenvalues
    # Normalize by total energy
    energy = eigenvalues / np.sum(eigenvalues)
    
    ax.semilogy(range(1, len(energy) + 1), energy, 'o-', linewidth = 1.5, markersize = 4)
    ax.set_xlabel('Mode index')
    ax.set_ylabel('Normalized energy')
    ax.set_title('POD Energy Spectrum')
    ax.grid(True, which = 'both', alpha = 0.3)
    
    fig.tight_layout()
    output_path = output_dir / 'plot_pod_spectrum.png'
    fig.savefig(output_path, dpi = DPI, bbox_inches = 'tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def plot_pod_modes(analyzer, x, y, output_dir):
    """Plot POD modes 1 and 2."""
    nx, ny = len(x), len(y)
    modes = analyzer.modes
    
    # Mode 1
    fig, ax = plt.subplots(figsize = (8, 6))
    mode1 = reshape_mode(modes[:, 0], nx, ny)
    plot_field(ax, x, y, mode1, 'POD Mode 1', cmap = 'RdBu_r')
    fig.tight_layout()
    output_path = output_dir / 'plot_pod_mode1.png'
    fig.savefig(output_path, dpi = DPI, bbox_inches = 'tight')
    plt.close(fig)
    print(f"Saved {output_path}")
    
    # Mode 2
    fig, ax = plt.subplots(figsize = (8, 6))
    mode2 = reshape_mode(modes[:, 1], nx, ny)
    plot_field(ax, x, y, mode2, 'POD Mode 2', cmap = 'RdBu_r')
    fig.tight_layout()
    output_path = output_dir / 'plot_pod_mode2.png'
    fig.savefig(output_path, dpi = DPI, bbox_inches = 'tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def plot_dmd_spectrum(analyzer, output_dir):
    """Plot DMD eigenvalues in complex plane with unit circle."""
    fig, ax = plt.subplots(figsize = (8, 8))
    
    eigenvalues = analyzer.eigenvalues
    
    # Unit circle
    theta = np.linspace(0, 2 * np.pi, 100)
    ax.plot(np.cos(theta), np.sin(theta), 'k--', linewidth = 1, alpha = 0.5, label = 'Unit circle')
    
    # Eigenvalues
    real_part = np.real(eigenvalues)
    imag_part = np.imag(eigenvalues)
    
    # Color by stability
    stable = np.abs(eigenvalues) <= 1.0
    ax.scatter(real_part[stable], imag_part[stable], c = 'blue', s = 50, 
               alpha = 0.6, label = 'Stable (|λ| ≤ 1)')
    ax.scatter(real_part[~stable], imag_part[~stable], c = 'red', s = 50, 
               alpha = 0.6, label = 'Unstable (|λ| > 1)')
    
    ax.set_xlabel('Re(λ)')
    ax.set_ylabel('Im(λ)')
    ax.set_title('DMD Eigenvalue Spectrum')
    ax.set_aspect('equal')
    ax.grid(True, alpha = 0.3)
    ax.legend(fontsize = 9)
    
    fig.tight_layout()
    output_path = output_dir / 'plot_dmd_spectrum.png'
    fig.savefig(output_path, dpi = DPI, bbox_inches = 'tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def plot_dmd_mode(analyzer, x, y, output_dir):
    """Plot leading DMD mode (real part)."""
    nx, ny = len(x), len(y)
    modes = analyzer.modes
    
    fig, ax = plt.subplots(figsize = (8, 6))
    mode1 = reshape_mode(np.real(modes[:, 0]), nx, ny)
    plot_field(ax, x, y, mode1, 'DMD Mode 1 (Real)', cmap = 'RdBu_r')
    fig.tight_layout()
    output_path = output_dir / 'plot_dmd_mode1.png'
    fig.savefig(output_path, dpi = DPI, bbox_inches = 'tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def plot_spod_spectrum(analyzer, output_dir):
    """Plot SPOD energy vs Strouhal number."""
    fig, ax = plt.subplots(figsize = (9, 6))
    
    St = analyzer.St
    eigenvalues = analyzer.eigenvalues  # shape (n_freq, n_modes)
    
    # Plot the first few eigenvalue curves
    n_curves = min(3, eigenvalues.shape[1])
    for i in range(n_curves):
        ax.semilogy(St, eigenvalues[:, i], 'o-', linewidth = 1.5, 
                   markersize = 3, label = f'Mode {i + 1}')
    
    ax.set_xlabel('Strouhal number (St)')
    ax.set_ylabel('SPOD energy')
    ax.set_title('SPOD Energy Spectrum vs Strouhal')
    ax.grid(True, which = 'both', alpha = 0.3)
    ax.legend(fontsize = 9)
    
    fig.tight_layout()
    output_path = output_dir / 'plot_spod_spectrum.png'
    fig.savefig(output_path, dpi = DPI, bbox_inches = 'tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def plot_spod_mode(analyzer, x, y, output_dir):
    """Plot leading SPOD mode at peak-energy frequency."""
    nx, ny = len(x), len(y)
    modes = analyzer.modes  # shape (n_freq, n_space, n_modes)
    eigenvalues = analyzer.eigenvalues  # shape (n_freq, n_modes)
    
    # Find frequency bin with highest energy (mode 0, all frequencies)
    peak_freq_idx = np.argmax(eigenvalues[:, 0])
    
    fig, ax = plt.subplots(figsize = (8, 6))
    mode1 = reshape_mode(np.real(modes[peak_freq_idx, :, 0]), nx, ny)
    st_peak = analyzer.St[peak_freq_idx]
    plot_field(ax, x, y, mode1, f'SPOD Mode 1 (St = {st_peak:.4f})', cmap = 'RdBu_r')
    fig.tight_layout()
    output_path = output_dir / 'plot_spod_mode1.png'
    fig.savefig(output_path, dpi = DPI, bbox_inches = 'tight')
    plt.close(fig)
    print(f"Saved {output_path}")


def main():
    """Main analysis pipeline."""
    print(f"Input NPZ: {NPZ_PATH}")
    print(f"Output dir: {CASE_DIR}")
    
    # Load grid
    x, y = load_grid()
    print(f"Grid: x ({len(x)},), y ({len(y)},)")
    
    # Use temporary directory for pyModal results (not for final PNGs)
    with tempfile.TemporaryDirectory() as temp_dir:
        # Run analyses
        pod = run_pod(str(NPZ_PATH), temp_dir)
        dmd = run_dmd(str(NPZ_PATH), temp_dir)
        spod = run_spod(str(NPZ_PATH), temp_dir)
        
        # Generate plots in case directory
        print("\nGenerating plots...")
        plot_pod_spectrum(pod, x, y, CASE_DIR)
        plot_pod_modes(pod, x, y, CASE_DIR)
        plot_dmd_spectrum(dmd, CASE_DIR)
        plot_dmd_mode(dmd, x, y, CASE_DIR)
        plot_spod_spectrum(spod, CASE_DIR)
        plot_spod_mode(spod, x, y, CASE_DIR)
    
    print("\nDone!")


if __name__ == '__main__':
    main()
