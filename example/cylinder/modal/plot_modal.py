#!/usr/bin/env python
"""
plot_modal.py: Unified plotting for modal analysis results (POD, DMD, SPOD)

Reads output files from nekStab modal_analysis.f90:
  - pod_energy.dat   : POD eigenvalue spectrum
  - dmd_spectrum.dat : DMD eigenvalues and frequencies
  - spod_spectrum.dat: SPOD energy vs frequency

Usage:
  python plot_modal.py              # auto-detect and plot all available
  python plot_modal.py --pod        # plot POD only
  python plot_modal.py --dmd        # plot DMD only
  python plot_modal.py --spod       # plot SPOD only

Output:
  pod_spectrum.png, dmd_spectrum.png, spod_spectrum.png
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
    """Read SPOD spectrum from spod_spectrum.dat"""

    def __init__(self, filename='spod_spectrum.dat'):
        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} not found')

        print(f'Reading {filename}')
        data = np.genfromtxt(filename, comments='#')

        if data.size == 0:
            raise ValueError(f'{filename} is empty')

        self.St = data[:, 0]                    # Strouhal numbers
        self.eigenvalues = data[:, 1:]          # eigenvalues (nfreq x nblk)
        self.nfreq = len(self.St)
        self.nblk = self.eigenvalues.shape[1]

        print(f'  {self.nfreq} frequencies, {self.nblk} SPOD modes per frequency')


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
    ax2.set_ylim(0, 105)

    # Reference lines at 90% and 99% with red markers on right axis
    for pct in [90, 99]:
        if cumulative[-1] >= pct:
            ax2.axhline(pct, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
            # Red marker on right edge (inside plot area)
            ax2.plot(n, pct, 's', color='firebrick', markersize=5, zorder=5)
            ax2.annotate(f'{pct}%', xy=(n - 0.5, pct), fontsize=7,
                         color='gray', va='center', ha='right')

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

    ax.set_xlabel(r'Strouhal number $St$')
    ax.set_ylabel('SPOD eigenvalue')
    ax.set_xlim(0, max(spod.St))
    ax.legend(loc='best', fontsize=7, ncol=2)
    ax.set_title('SPOD Spectrum')
    ax.grid(True, which='major', linestyle='-', linewidth=0.3, alpha=0.5)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.2, alpha=0.3)

    fig.tight_layout()

    fname = f'spod_spectrum.{FORMAT}'
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

    plots_made = 0

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

    # SPOD
    if plot_spod_flag:
        try:
            spod = SPODSpectrum()
            plot_spod(spod)
            plots_made += 1
        except (FileNotFoundError, ValueError) as e:
            print(f'  Skipping SPOD: {e}')

    print('-' * 50)
    if plots_made == 0:
        print('No data files found. Run modal analysis first.')
        print('Expected files: pod_energy.dat, dmd_spectrum.dat, spod_spectrum.dat')
    else:
        print(f'Done. Generated {plots_made} plot(s).')
