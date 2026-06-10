#!/usr/bin/env python3
"""plot.py: Visualize Poiseuille OTD results.

OUTPUTS: plot_lyapunov_exponents.png, plot_otd_residuals.png, plot_otd_mode1.png, plot_otd_mode2.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def main():
    nk.configure_style()

    has_dat = any((CASE_DIR / f).exists()
                  for f in ('otd_ftle.dat', 'otd_eigenvalues.dat', 'otd_growth_rates.dat'))
    has_resid = (CASE_DIR / 'otd_residuals.dat').exists()
    ff1 = nk.find_fields('r01poiseuille_OTD0.f00001', CASE_DIR)
    ff2 = nk.find_fields('r02poiseuille_OTD0.f00001', CASE_DIR)

    if has_dat:
        fig, ax_otd = plt.subplots(
            1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_otd_exponents(ax_otd, CASE_DIR)
        ax_otd.set_title('Lyapunov exponents', fontsize=8)
        output = CASE_DIR / 'plot_lyapunov_exponents.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    if has_resid:
        fig, ax_res = plt.subplots(
            1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        data = np.genfromtxt(str(CASE_DIR / 'otd_residuals.dat'))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        t = data[:, 0]
        n_modes = data.shape[1] - 1
        colors_cycle = ['b', 'g', 'r', 'c', 'm']
        for j in range(n_modes):
            c = colors_cycle[j % len(colors_cycle)]
            vals = data[:, j + 1]
            mask = vals > 0  # skip zeros for log scale
            ax_res.semilogy(t[mask], vals[mask], lw=0.8, color=c,
                            label=rf'$|\Delta\lambda_{j+1}|$')
        ax_res.axhline(1e-6, color='k', ls='--', lw=0.6,
                       label=r'tol $= 10^{-6}$')
        ax_res.set_xlabel('$t$')
        ax_res.set_ylabel(r'$|\Delta \mathrm{FTLE}|$')
        ax_res.legend(fontsize=5, loc='best')
        ax_res.grid(True, ls=':', lw=0.3, alpha=0.5)
        output = CASE_DIR / 'plot_otd_residuals.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    for ff, title, filename in [
        (ff1, 'OTD mode 1 ($v_y$)', 'plot_otd_mode1.png'),
        (ff2, 'OTD mode 2 ($v_y$)', 'plot_otd_mode2.png'),
    ]:
        if ff:
            fig, ax_m = plt.subplots(
                1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
            x, y, fields, _ = nk.read_field(ff[0])
            triang = nk.make_triangulation(x, y)
            q = fields.get('vy', fields.get('vx'))
            if q is not None:
                bd = np.nanpercentile(np.abs(q), 99)
                if bd < 1e-15:
                    bd = 1.0
                cf = nk.tricontourf(ax_m, triang, q, levels=257,
                                    cmap='RdBu', vmin=-bd, vmax=bd,
                                    extend='both')
                bdr = round(bd, 2)
                nk.inset_colorbar(ax_m, cf, orientation='horizontal',
                                  width="50%", height="5%", loc=9,
                                  ticks=[-bdr, bdr],
                                  tick_labels=[f'{-bdr}', f'{bdr}'])
            ax_m.set_aspect('equal')
            ax_m.set_xlabel(r'$x$', labelpad=-1)
            ax_m.set_ylabel(r'$y$', labelpad=1)
            ax_m.spines['right'].set_visible(False)
            ax_m.spines['top'].set_visible(False)
            ax_m.set_title(title, fontsize=8)
            output = CASE_DIR / filename
            fig.savefig(output, dpi=600, bbox_inches='tight')
            print(f'Saved {output}')
            plt.close(fig)


if __name__ == '__main__':
    main()
