#!/usr/bin/env python3
"""plot.py: Visualize OTD analysis results for cylinder Re=180.

Panels (emitted only when data is present):
  plot_lyapunov_exponents.png - Lyapunov exponents / FTLE
  plot_otd_residuals.png      - OTD convergence residuals
  plot_otd_mode1.png          - OTD mode 1 (v_y)
  plot_otd_mode2.png          - OTD mode 2 (v_y)

OUTPUTS: plot_lyapunov_exponents.png, plot_otd_residuals.png, plot_otd_mode1.png, plot_otd_mode2.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def find_otd_mode(mode_index):
    """Find OTD mode field file for cylinder. Tries multiple naming patterns."""
    patterns = [
        f'ip{mode_index}*1cyl0.f00001',
        f'r0{mode_index + 1}*1cyl0.f00001',
        f'r{mode_index + 1:02d}*1cyl0.f00001',
    ]
    for pat in patterns:
        matches = sorted(CASE_DIR.glob(pat))
        if matches:
            return matches[0]
    return None


def main():
    nk.configure_style()

    # PANEL: LYAPUNOV EXPONENTS
    has_dat = any((CASE_DIR / f).exists()
                  for f in ('otd_ftle.dat', 'otd_eigenvalues.dat', 'otd_growth_rates.dat'))
    if has_dat:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_otd_exponents(ax, CASE_DIR)
        ax.set_title('Lyapunov exponents', fontsize=8)
        output = CASE_DIR / 'plot_lyapunov_exponents.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # PANEL: OTD RESIDUALS
    resid_file = CASE_DIR / 'otd_residuals.dat'
    if resid_file.exists():
        data = np.genfromtxt(str(resid_file))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        t = data[:, 0]
        n_modes = data.shape[1] - 1
        colors_cycle = ['b', 'g', 'r', 'c', 'm']
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        for j in range(n_modes):
            c = colors_cycle[j % len(colors_cycle)]
            vals = data[:, j + 1]
            mask = vals > 0  # skip zeros for log scale
            ax.semilogy(t[mask], vals[mask], lw=0.8, color=c,
                        label=rf'$|\Delta\lambda_{j+1}|$')
        ax.axhline(1e-6, color='k', ls='--', lw=0.6,
                   label=r'tol $= 10^{-6}$')
        ax.set_xlabel('$t$')
        ax.set_ylabel(r'$|\Delta \mathrm{FTLE}|$')
        ax.legend(fontsize=5, loc='best')
        ax.grid(True, ls=':', lw=0.3, alpha=0.5)
        output = CASE_DIR / 'plot_otd_residuals.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # PANEL: OTD MODE 1 and MODE 2
    for mode_idx, title, fname in [
        (0, 'OTD mode 1 ($v_y$)', 'plot_otd_mode1.png'),
        (1, 'OTD mode 2 ($v_y$)', 'plot_otd_mode2.png'),
    ]:
        ff = find_otd_mode(mode_idx)
        if ff is None:
            continue
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        x, y, fields, _ = nk.read_field(ff)
        triang = nk.make_triangulation(x, y)
        q = fields.get('vy', fields.get('vx'))
        if q is not None:
            bd = np.nanpercentile(np.abs(q), 99)
            if bd < 1e-15:
                bd = 1.0
            cf = nk.tricontourf(ax, triang, q, levels=257,
                                cmap='RdBu', vmin=-bd, vmax=bd, extend='both')
            bdr = round(bd, 2)
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[-bdr, bdr],
                              tick_labels=[f'{-bdr}', f'{bdr}'])
        ax.set_aspect('equal')
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-2, 20)
        ax.set_ylim(-4, 4)
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(title, fontsize=8)
        output = CASE_DIR / fname
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)


if __name__ == '__main__':
    main()
