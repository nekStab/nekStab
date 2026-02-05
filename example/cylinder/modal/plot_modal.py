#!/usr/bin/env python
"""
plot_modal.py: Unified plotting for modal analysis results (POD, DMD, SPOD)

Reads output files from nekStab modal_analysis.f90:
  - pod_energy.dat          : POD eigenvalue spectrum
  - pod_fft_spectrum.dat    : POD-FFT power spectra
  - dmd_spectrum.dat        : DMD eigenvalues and frequencies
  - spod_spectrum.dat       : Batch SPOD spectrum (legacy)
  - spod_stream_spectrum.dat: Streaming SPOD spectrum

Usage:
  python plot_modal.py              # auto-detect and plot all available
  python plot_modal.py --pod        # plot POD only
  python plot_modal.py --dmd        # plot DMD only
  python plot_modal.py --spod       # plot SPOD only
  python plot_modal.py --compare    # plot POD-FFT vs SPOD comparison

Output:
  pod_spectrum.png, dmd_spectrum.png, spod_spectrum.png,
  pod_fft_spectrum.png, spectral_comparison.png
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# -----------------------------------------------------------------------------
# Matplotlib configuration (consistent with nekStab style)
# -----------------------------------------------------------------------------
params = {
    'text.usetex': False,
    'font.size': 9,
    'legend.fontsize': 8,
    'legend.handlelength': 1.5,
    'axes.labelsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
}
plt.rcParams.update(params)

FORMAT = 'png'
DPI = 600
BBOX = 'tight'

# Reference Strouhal number for Re=100 cylinder (from stability analysis)
ST_REF = 0.16594721970048534


# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------

class PODSpectrum:
    """Read POD eigenvalue spectrum from pod_energy.dat"""

    def __init__(self, filename='pod_energy.dat'):
        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} not found')

        print(f'Reading {filename}')
        data = np.genfromtxt(filename, comments='#')

        self.mode = data[:, 0].astype(int)
        self.eigenvalue = data[:, 1]
        self.energy_pct = data[:, 2]
        self.cumulative_pct = data[:, 3]
        self.nmodes = len(self.mode)

        print(f'  {self.nmodes} modes, total energy captured: {self.cumulative_pct[-1]:.2f}%')


class DMDSpectrum:
    """Read DMD eigenvalue spectrum from dmd_spectrum.dat"""

    def __init__(self, filename='dmd_spectrum.dat'):
        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} not found')

        print(f'Reading {filename}')
        data = np.genfromtxt(filename, comments='#')

        self.mode = data[:, 0].astype(int)
        self.mu_abs = data[:, 1]      # |mu|
        self.sigma = data[:, 2]       # growth rate
        self.omega = data[:, 3]       # angular frequency
        self.St = data[:, 4]          # Strouhal number
        self.mu_real = data[:, 5]     # Re(mu)
        self.mu_imag = data[:, 6]     # Im(mu)
        self.nmodes = len(self.mode)

        # Find unstable modes
        n_unstable = np.sum(self.sigma > 0)
        print(f'  {self.nmodes} modes, {n_unstable} unstable (sigma > 0)')


class SPODSpectrum:
    """Read SPOD spectrum from spod_spectrum.dat or spod_stream_spectrum.dat"""

    def __init__(self, filename=None):
        # Auto-detect file: prefer streaming SPOD
        if filename is None:
            if os.path.exists('spod_stream_spectrum.dat'):
                filename = 'spod_stream_spectrum.dat'
                self.method = 'streaming'
            elif os.path.exists('spod_spectrum.dat'):
                filename = 'spod_spectrum.dat'
                self.method = 'batch'
            else:
                raise FileNotFoundError('No SPOD spectrum file found')
        else:
            self.method = 'streaming' if 'stream' in filename else 'batch'

        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} not found')

        print(f'Reading {filename} ({self.method} SPOD)')
        data = np.genfromtxt(filename, comments='#')

        if data.size == 0:
            raise ValueError(f'{filename} is empty')

        self.filename = filename
        self.St = data[:, 0]                    # Strouhal numbers
        self.eigenvalues = data[:, 1:]          # eigenvalues (nfreq x nblk)
        self.nfreq = len(self.St)
        self.nblk = self.eigenvalues.shape[1]

        # Find peak frequency
        idx = np.argmax(self.eigenvalues[1:, 0]) + 1  # Skip DC
        self.peak_St = self.St[idx]
        self.peak_lambda = self.eigenvalues[idx, 0]

        print(f'  {self.nfreq} frequencies, {self.nblk} SPOD modes per frequency')
        print(f'  Peak: St = {self.peak_St:.4f}, λ₁ = {self.peak_lambda:.4f}')


class PODFFTSpectrum:
    """Read POD-FFT spectrum from pod_fft_spectrum.dat"""

    def __init__(self, filename='pod_fft_spectrum.dat'):
        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} not found')

        print(f'Reading {filename}')
        data = np.genfromtxt(filename, comments='#')

        if data.size == 0:
            raise ValueError(f'{filename} is empty')

        self.St = data[:, 0]                    # Strouhal numbers
        self.power = data[:, 1:]                # power spectra (nfreq x nmodes)
        self.nfreq = len(self.St)
        self.nmodes = self.power.shape[1]

        # Find peak frequency for each mode (excluding DC)
        self.peak_St = []
        self.peak_power = []
        for i in range(self.nmodes):
            idx = np.argmax(self.power[1:, i]) + 1  # Skip DC (index 0)
            self.peak_St.append(self.St[idx])
            self.peak_power.append(self.power[idx, i])

        print(f'  {self.nfreq} frequencies, {self.nmodes} POD modes')
        print(f'  Peak St for mode 1: {self.peak_St[0]:.4f}')


# -----------------------------------------------------------------------------
# Plotting functions
# -----------------------------------------------------------------------------

def plot_pod(pod, max_modes=20):
    """Plot POD energy spectrum: bar chart + cumulative line"""

    n = min(max_modes, pod.nmodes)
    modes = pod.mode[:n]
    energy = pod.energy_pct[:n]
    cumulative = pod.cumulative_pct[:n]

    fig, ax1 = plt.subplots(figsize=(5, 3.5))

    # Bar chart for individual mode energy
    ax1.bar(modes, energy, color='steelblue', edgecolor='k',
            linewidth=0.5, alpha=0.8, label='Mode energy')
    ax1.set_xlabel('POD mode')
    ax1.set_ylabel('Energy (%)', color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.set_xlim(0.5, n + 0.5)
    ax1.set_ylim(0, max(energy) * 1.15)

    # X-axis: integer ticks only
    ax1.set_xticks(modes[::2])  # every other mode to avoid crowding
    ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))

    # Cumulative energy on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(modes, cumulative, 'o-', color='firebrick', markersize=4,
             linewidth=1.2, label='Cumulative')
    ax2.set_ylabel('Cumulative energy (%)', color='firebrick')
    ax2.tick_params(axis='y', labelcolor='firebrick')
    ax2.set_ylim(0, 100)

    # Reference lines at 99% and 90% (order matters for legend)
    if cumulative[-1] >= 99:
        ax2.axhline(99, color='green', linestyle='--', linewidth=1.0, label='99%')
    if cumulative[-1] >= 90:
        ax2.axhline(90, color='orange', linestyle='--', linewidth=1.0, label='90%')

    ax2.legend(loc='lower right', fontsize=7)
    ax1.set_title('POD Energy Spectrum')
    fig.tight_layout()

    fname = f'pod_spectrum.{FORMAT}'
    plt.savefig(fname, format=FORMAT, dpi=DPI, bbox_inches=BBOX)
    print(f'Saved {fname}')
    plt.close()


def plot_dmd(dmd):
    """Plot DMD spectrum: complex plane + growth rate vs frequency"""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))

    # --- Left panel: Complex plane (mu) ---
    theta = np.linspace(0, 2*np.pi, 200)
    ax1.plot(np.cos(theta), np.sin(theta), 'r-', linewidth=0.8,
             label='Unit circle')

    # Color by stability: blue=stable, red=unstable
    stable = dmd.mu_abs <= 1.0
    unstable = dmd.mu_abs > 1.0

    ax1.scatter(dmd.mu_real[stable], dmd.mu_imag[stable],
                s=25, c='steelblue', edgecolors='k', linewidth=0.3,
                alpha=0.8, label=f'Stable ({np.sum(stable)})', zorder=3)

    if np.any(unstable):
        ax1.scatter(dmd.mu_real[unstable], dmd.mu_imag[unstable],
                    s=40, c='firebrick', edgecolors='k', linewidth=0.3,
                    marker='D', alpha=0.9, label=f'Unstable ({np.sum(unstable)})',
                    zorder=4)

    ax1.axhline(0, color='k', linewidth=0.3, linestyle=':')
    ax1.axvline(0, color='k', linewidth=0.3, linestyle=':')
    ax1.set_xlabel(r'Re$(\mu)$')
    ax1.set_ylabel(r'Im$(\mu)$')
    ax1.set_aspect('equal')
    ax1.set_xlim(-1.3, 1.3)
    ax1.set_ylim(-1.3, 1.3)
    ax1.legend(loc='upper left', fontsize=7)
    ax1.set_title('DMD Eigenvalues')

    # --- Right panel: Growth rate vs Strouhal ---
    ax2.scatter(np.abs(dmd.St), dmd.sigma, s=25, c='steelblue',
                edgecolors='k', linewidth=0.3, alpha=0.8)

    # Highlight unstable modes
    if np.any(dmd.sigma > 0):
        ax2.scatter(np.abs(dmd.St[dmd.sigma > 0]), dmd.sigma[dmd.sigma > 0],
                    s=40, c='firebrick', edgecolors='k', linewidth=0.3,
                    marker='D', alpha=0.9, zorder=4)

    ax2.axhline(0, color='r', linewidth=0.8, linestyle='-', label=r'$\sigma=0$')

    # Reference Strouhal
    ax2.axvline(ST_REF, color='green', linestyle='--', linewidth=1.0,
                alpha=0.8, label=rf'$St_{{DNS}} = {ST_REF:.3f}$')

    ax2.set_xlabel(r'Strouhal number $St$')
    ax2.set_ylabel(r'Growth rate $\sigma$')
    ax2.set_xlim(left=0)
    ax2.legend(loc='best', fontsize=7)
    ax2.set_title('DMD Growth Rates')

    fig.tight_layout()

    fname = f'dmd_spectrum.{FORMAT}'
    plt.savefig(fname, format=FORMAT, dpi=DPI, bbox_inches=BBOX)
    print(f'Saved {fname}')
    plt.close()


def plot_spod(spod, max_modes=5):
    """Plot SPOD spectrum: energy vs Strouhal for each mode rank"""

    fig, ax = plt.subplots(figsize=(5.5, 3.5))

    # Color palette for different mode ranks
    colors = plt.cm.viridis(np.linspace(0, 0.85, min(max_modes, spod.nblk)))

    n_to_plot = min(max_modes, spod.nblk)
    for i in range(n_to_plot):
        eig = spod.eigenvalues[:, i]
        # Avoid log(0) issues
        eig = np.maximum(eig, 1e-20)

        lw = 1.5 if i == 0 else 0.8
        alpha = 1.0 if i == 0 else 0.7
        ax.semilogy(spod.St, eig, '-', color=colors[i], linewidth=lw,
                    alpha=alpha, label=f'Mode {i+1}')

    # Reference line at peak
    ax.axvline(spod.peak_St, color='red', linestyle='--', linewidth=1.0,
               alpha=0.7, label=f'Peak St = {spod.peak_St:.3f}')

    # Reference Strouhal for Re=100 cylinder
    ax.axvline(ST_REF, color='green', linestyle=':', linewidth=1.0,
               alpha=0.7, label=f'DNS St = {ST_REF:.3f}')

    ax.set_xlabel(r'Strouhal number $St$')
    ax.set_ylabel('SPOD eigenvalue')
    ax.set_xlim(0, max(spod.St))
    ax.legend(loc='best', fontsize=7, ncol=2)
    title = f'SPOD Spectrum ({spod.method.capitalize()})'
    ax.set_title(title)
    ax.grid(True, which='major', linestyle='-', linewidth=0.3, alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.2, alpha=0.3)

    fig.tight_layout()

    fname = f'spod_spectrum.{FORMAT}'
    plt.savefig(fname, format=FORMAT, dpi=DPI, bbox_inches=BBOX)
    print(f'Saved {fname}')
    plt.close()


def plot_pod_fft(pod_fft, max_modes=5):
    """Plot POD-FFT spectrum: power vs Strouhal for each POD mode"""

    fig, ax = plt.subplots(figsize=(5.5, 3.5))

    # Color palette
    colors = plt.cm.tab10(np.linspace(0, 1, min(max_modes, pod_fft.nmodes)))

    n_to_plot = min(max_modes, pod_fft.nmodes)
    for i in range(n_to_plot):
        power = pod_fft.power[:, i]
        # Avoid log(0) issues
        power = np.maximum(power, 1e-20)

        lw = 1.8 if i < 2 else 1.0
        alpha = 1.0 if i < 2 else 0.7
        ax.semilogy(pod_fft.St, power, 'o-', color=colors[i], linewidth=lw,
                    markersize=4, alpha=alpha, label=f'POD mode {i+1}')

    # Reference Strouhal
    ax.axvline(ST_REF, color='red', linestyle='--', linewidth=1.0,
               alpha=0.7, label=rf'$St = {ST_REF:.4f}$ (Re=100)')
    # Add marker at the peak of mode 1 if it's near the reference
    if pod_fft.nmodes > 0:
        y_pos = pod_fft.peak_power[0]
        ax.plot(ST_REF, y_pos, 'r*', markersize=12, markeredgecolor='k',
                markeredgewidth=0.5, zorder=5)

    ax.set_xlabel(r'Strouhal number $St$')
    ax.set_ylabel('Power spectral density')
    ax.set_xlim(0, min(1.0, max(pod_fft.St)))
    ax.legend(loc='best', fontsize=7)
    ax.set_title('POD-FFT Spectral Analysis')
    ax.grid(True, which='major', linestyle='-', linewidth=0.3, alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.2, alpha=0.3)

    fig.tight_layout()

    fname = f'pod_fft_spectrum.{FORMAT}'
    plt.savefig(fname, format=FORMAT, dpi=DPI, bbox_inches=BBOX)
    print(f'Saved {fname}')
    plt.close()


def plot_spectral_comparison(pod_fft=None, spod_stream=None, spod_batch=None):
    """Plot comparison of POD-FFT, Batch SPOD, and Streaming SPOD spectra"""

    if pod_fft is None and spod_stream is None and spod_batch is None:
        print('  No data for comparison plot')
        return

    fig, ax = plt.subplots(figsize=(6, 4))

    # POD-FFT (modes 1-2)
    if pod_fft is not None:
        power1 = np.maximum(pod_fft.power[:, 0], 1e-20)
        ax.semilogy(pod_fft.St, power1, 'b-', linewidth=1.5,
                    label='POD-FFT mode 1', alpha=0.9)
        if pod_fft.nmodes > 1:
            power2 = np.maximum(pod_fft.power[:, 1], 1e-20)
            ax.semilogy(pod_fft.St, power2, 'b--', linewidth=1.2,
                        label='POD-FFT mode 2', alpha=0.7)

    # Batch SPOD (if available)
    if spod_batch is not None:
        eig1 = np.maximum(spod_batch.eigenvalues[:, 0], 1e-20)
        ax.semilogy(spod_batch.St, eig1, 'g-', linewidth=1.5,
                    label='Batch SPOD mode 1', alpha=0.9)

    # Streaming SPOD
    if spod_stream is not None:
        eig1 = np.maximum(spod_stream.eigenvalues[:, 0], 1e-20)
        ax.semilogy(spod_stream.St, eig1, 'r-', linewidth=1.5,
                    label='Streaming SPOD mode 1', alpha=0.9)
        if spod_stream.nblk > 1:
            eig2 = np.maximum(spod_stream.eigenvalues[:, 1], 1e-20)
            ax.semilogy(spod_stream.St, eig2, 'r--', linewidth=1.2,
                        label='Streaming SPOD mode 2', alpha=0.7)

    # Reference lines
    ax.axvline(ST_REF, color='black', linestyle=':', linewidth=1.5,
               alpha=0.8, label=f'DNS St = {ST_REF:.3f}')

    ax.set_xlabel(r'Strouhal number $St$')
    ax.set_ylabel('Spectral energy / Power')
    ax.set_xlim(0, 0.8)
    ax.legend(loc='upper right', fontsize=7)
    ax.set_title('Spectral Analysis Comparison')
    ax.grid(True, which='major', linestyle='-', linewidth=0.3, alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.2, alpha=0.3)

    # Add annotation for peak
    spod = spod_stream or spod_batch
    if spod is not None:
        ax.annotate(f'Peak: St={spod.peak_St:.3f}',
                    xy=(spod.peak_St, spod.peak_lambda),
                    xytext=(spod.peak_St + 0.1, spod.peak_lambda),
                    fontsize=8, color='red',
                    arrowprops=dict(arrowstyle='->', color='red', lw=0.8))

    fig.tight_layout()

    fname = f'spectral_comparison.{FORMAT}'
    plt.savefig(fname, format=FORMAT, dpi=DPI, bbox_inches=BBOX)
    print(f'Saved {fname}')
    plt.close()


def plot_spod_validation(spod_batch, spod_stream):
    """Plot validation: Batch SPOD vs Streaming SPOD

    This comparison validates that the streaming algorithm gives the same
    results as the (slower) reference batch implementation.
    """

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # --- Left panel: Eigenvalue comparison ---
    eig_batch = np.maximum(spod_batch.eigenvalues[:, 0], 1e-20)
    eig_stream = np.maximum(spod_stream.eigenvalues[:, 0], 1e-20)

    ax1.semilogy(spod_batch.St, eig_batch, 'g-o', linewidth=1.5,
                 markersize=5, label='Batch SPOD', alpha=0.9)
    ax1.semilogy(spod_stream.St, eig_stream, 'r--s', linewidth=1.5,
                 markersize=4, label='Streaming SPOD', alpha=0.8)

    ax1.axvline(ST_REF, color='black', linestyle=':', linewidth=1.0,
                alpha=0.6, label=f'DNS St = {ST_REF:.3f}')

    ax1.set_xlabel(r'Strouhal number $St$')
    ax1.set_ylabel('SPOD eigenvalue (mode 1)')
    ax1.legend(loc='best', fontsize=8)
    ax1.set_title('Eigenvalue Comparison')
    ax1.grid(True, which='major', linestyle='-', linewidth=0.3, alpha=0.5)

    # --- Right panel: Relative error ---
    # Interpolate to common grid if needed
    if len(spod_batch.St) == len(spod_stream.St):
        rel_error = np.abs(eig_stream - eig_batch) / np.maximum(eig_batch, 1e-20)
        ax2.semilogy(spod_batch.St, rel_error * 100, 'k-o', linewidth=1.2,
                     markersize=4)
        ax2.set_xlabel(r'Strouhal number $St$')
        ax2.set_ylabel('Relative error (%)')
        ax2.set_title('Streaming vs Batch Error')
        ax2.grid(True, which='major', linestyle='-', linewidth=0.3, alpha=0.5)

        max_err = np.max(rel_error) * 100
        ax2.axhline(max_err, color='red', linestyle='--', linewidth=0.8,
                    label=f'Max error: {max_err:.2f}%')
        ax2.legend(loc='best', fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'Different frequency grids\nCannot compute error',
                 ha='center', va='center', transform=ax2.transAxes)

    fig.tight_layout()

    fname = f'spod_validation.{FORMAT}'
    plt.savefig(fname, format=FORMAT, dpi=DPI, bbox_inches=BBOX)
    print(f'Saved {fname}')
    plt.close()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    print('=' * 50)
    print('  Modal Analysis Spectrum Plotter')
    print('=' * 50)

    # Parse arguments
    plot_all = len(sys.argv) == 1
    plot_pod_flag = plot_all or '--pod' in sys.argv
    plot_dmd_flag = plot_all or '--dmd' in sys.argv
    plot_spod_flag = plot_all or '--spod' in sys.argv
    plot_compare_flag = plot_all or '--compare' in sys.argv

    plots_made = 0
    pod_fft_data = None
    spod_data = None

    # POD
    if plot_pod_flag:
        try:
            pod = PODSpectrum()
            plot_pod(pod)
            plots_made += 1
        except FileNotFoundError as e:
            print(f'  Skipping POD: {e}')

    # DMD
    if plot_dmd_flag:
        try:
            dmd = DMDSpectrum()
            plot_dmd(dmd)
            plots_made += 1
        except FileNotFoundError as e:
            print(f'  Skipping DMD: {e}')

    # SPOD - try both streaming and batch
    spod_stream = None
    spod_batch = None

    if plot_spod_flag:
        # Try streaming SPOD first (preferred)
        try:
            spod_stream = SPODSpectrum('spod_stream_spectrum.dat')
            plot_spod(spod_stream)
            plots_made += 1
        except (FileNotFoundError, ValueError) as e:
            print(f'  Skipping streaming SPOD: {e}')

        # Try batch SPOD
        try:
            spod_batch = SPODSpectrum('spod_spectrum.dat')
            if spod_stream is None:
                plot_spod(spod_batch)
                plots_made += 1
        except (FileNotFoundError, ValueError) as e:
            print(f'  Skipping batch SPOD: {e}')

    # POD-FFT
    try:
        pod_fft_data = PODFFTSpectrum()
        plot_pod_fft(pod_fft_data)
        plots_made += 1
    except (FileNotFoundError, ValueError) as e:
        print(f'  Skipping POD-FFT: {e}')

    # Comparison plot: POD-FFT vs SPOD(s)
    if plot_compare_flag:
        if pod_fft_data is not None or spod_stream is not None or spod_batch is not None:
            plot_spectral_comparison(pod_fft_data, spod_stream, spod_batch)
            plots_made += 1

    # Validation plot: Batch SPOD vs Streaming SPOD (if both exist)
    if spod_batch is not None and spod_stream is not None:
        print('  Found both batch and streaming SPOD - creating validation plot')
        plot_spod_validation(spod_batch, spod_stream)
        plots_made += 1

    print('-' * 50)
    if plots_made == 0:
        print('No data files found. Run modal analysis first.')
        print('Expected files: pod_energy.dat, dmd_spectrum.dat,')
        print('                spod_stream_spectrum.dat, pod_fft_spectrum.dat')
    else:
        print(f'Done. Generated {plots_made} plot(s).')
