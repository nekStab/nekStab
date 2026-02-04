#!/usr/bin/env python3
"""
Plot modal analysis results for cylinder Re=100.

Expected: Vortex shedding at St ≈ 0.164-0.167
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Expected Strouhal number range for Re=100
ST_MIN, ST_MAX = 0.160, 0.170


def plot_pod_spectrum():
    """Plot POD eigenvalue spectrum."""
    try:
        data = np.loadtxt('pod_energy.dat', comments='#')
    except FileNotFoundError:
        print("pod_energy.dat not found")
        return

    mode = data[:, 0]
    eigval = data[:, 1]
    cumsum = data[:, 3]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    # Eigenvalue spectrum
    ax1.semilogy(mode, eigval, 'bo-', markersize=4)
    ax1.set_xlabel('Mode')
    ax1.set_ylabel('Eigenvalue')
    ax1.set_title('POD Eigenvalue Spectrum')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, min(50, len(mode)))

    # Cumulative energy
    ax2.plot(mode, cumsum, 'r-', linewidth=2)
    ax2.axhline(90, color='k', linestyle='--', alpha=0.5, label='90%')
    ax2.axhline(99, color='k', linestyle=':', alpha=0.5, label='99%')
    ax2.set_xlabel('Mode')
    ax2.set_ylabel('Cumulative Energy (%)')
    ax2.set_title('POD Cumulative Energy')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0, min(50, len(mode)))
    ax2.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig('pod_spectrum.png', dpi=150)
    print("Saved pod_spectrum.png")
    plt.close()


def plot_dmd_spectrum():
    """Plot DMD eigenvalue spectrum."""
    try:
        data = np.loadtxt('dmd_spectrum.dat', comments='#')
    except FileNotFoundError:
        print("dmd_spectrum.dat not found")
        return

    mu_mag = data[:, 1]
    sigma = data[:, 2]  # Growth rate
    St = data[:, 4]     # Strouhal number (column header is St)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    # Eigenvalue magnitude vs Strouhal number
    ax1.scatter(St, mu_mag, c='b', s=30, alpha=0.7)
    ax1.axvline(ST_MIN, color='r', linestyle='--', alpha=0.5)
    ax1.axvline(ST_MAX, color='r', linestyle='--', alpha=0.5,
                label=f'Expected St={ST_MIN:.3f}-{ST_MAX:.3f}')
    ax1.axhline(1.0, color='k', linestyle='-', alpha=0.3)
    ax1.set_xlabel('St')
    ax1.set_ylabel('|μ|')
    ax1.set_title('DMD Eigenvalue Magnitude')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Growth rate vs Strouhal number
    ax2.scatter(St, sigma, c='b', s=30, alpha=0.7)
    ax2.axvline(ST_MIN, color='r', linestyle='--', alpha=0.5)
    ax2.axvline(ST_MAX, color='r', linestyle='--', alpha=0.5)
    ax2.axhline(0.0, color='k', linestyle='-', alpha=0.3)
    ax2.set_xlabel('St')
    ax2.set_ylabel('Growth rate σ')
    ax2.set_title('DMD Growth Rate')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('dmd_spectrum.png', dpi=150)
    print("Saved dmd_spectrum.png")
    plt.close()

    # Find dominant mode
    idx = np.argmax(mu_mag)
    print(f"DMD: Dominant mode at St = {St[idx]:.4f}, |μ| = {mu_mag[idx]:.4f}")


def plot_spod_spectrum():
    """Plot SPOD spectrum."""
    try:
        data = np.loadtxt('spod_spectrum.dat', comments='#')
    except FileNotFoundError:
        print("spod_spectrum.dat not found")
        return

    St = data[:, 0]  # Strouhal number (column header is St)
    # Eigenvalues are in columns 1, 2, 3, ...
    evals = data[:, 1:]

    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot first few eigenvalues
    nplot = min(5, evals.shape[1])
    for i in range(nplot):
        ax.semilogy(St, evals[:, i], label=f'Mode {i+1}', alpha=0.8)

    # Mark expected shedding frequency
    ax.axvline(ST_MIN, color='r', linestyle='--', alpha=0.5)
    ax.axvline(ST_MAX, color='r', linestyle='--', alpha=0.5,
               label=f'Expected St={ST_MIN:.3f}-{ST_MAX:.3f}')

    ax.set_xlabel('St')
    ax.set_ylabel('SPOD Eigenvalue')
    ax.set_title('SPOD Spectrum')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('spod_spectrum.png', dpi=150)
    print("Saved spod_spectrum.png")
    plt.close()

    # Find peak frequency
    idx = np.argmax(evals[:, 0])
    print(f"SPOD: Peak at St = {St[idx]:.4f}")


def main():
    print("Cylinder Re=100 Modal Analysis Results")
    print("=" * 50)
    print(f"Expected: Vortex shedding at St ≈ {ST_MIN}-{ST_MAX}")
    print("")

    plot_pod_spectrum()
    plot_dmd_spectrum()
    plot_spod_spectrum()

    print("")
    print("Done! Check generated PNG files.")


if __name__ == '__main__':
    main()
