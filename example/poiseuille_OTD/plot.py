#!/usr/bin/env python3
"""plot.py: Visualize Poiseuille OTD results.

OUTPUTS: plot.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OUTPUT = CASE_DIR / 'plot.png'


def main():
    nk.configure_style()

    has_dat = any((CASE_DIR / f).exists()
                  for f in ('otd_ftle.dat', 'otd_eigenvalues.dat', 'otd_growth_rates.dat'))
    has_resid = (CASE_DIR / 'otd_residuals.dat').exists()
    ff1 = nk.find_fields('r01poiseuille_OTD0.f00001', CASE_DIR)
    ff2 = nk.find_fields('r02poiseuille_OTD0.f00001', CASE_DIR)

    # Determine layout
    n_left = int(has_dat) + int(has_resid)
    has_mode = bool(ff1)

    if not has_dat and not has_resid and not has_mode:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        ax.text(0.5, 0.5, 'No OTD data found',
                transform=ax.transAxes, ha='center')
        fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
        plt.close()
        return

    if n_left == 2 and has_mode:
        # 4 panels: exponents + convergence left, 2 modes right
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.9))
        gs = fig.add_gridspec(2, 2, width_ratios=[0.8, 1],
                              hspace=0.15, wspace=0.35)
        ax_otd = fig.add_subplot(gs[0, 0])
        ax_res = fig.add_subplot(gs[1, 0])
        ax_mode1 = fig.add_subplot(gs[0, 1])
        ax_mode2 = fig.add_subplot(gs[1, 1])
    elif n_left == 2:
        fig, (ax_otd, ax_res) = plt.subplots(
            1, 2, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.5))
        ax_mode = None
    elif has_dat and has_mode:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_otd = fig.add_subplot(gs[0])
        ax_mode1 = fig.add_subplot(gs[1])
        ax_mode2 = None
        ax_res = None
    elif has_dat:
        fig, ax_otd = plt.subplots(
            1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        ax_res = None
        ax_mode1 = None
        ax_mode2 = None
    else:
        fig, ax_res = plt.subplots(
            1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        ax_otd = None
        ax_mode1 = None
        ax_mode2 = None

    labels = iter('abcdefgh')

    if has_dat and ax_otd is not None:
        nk.plot_otd_exponents(ax_otd, CASE_DIR)
        ax_otd.set_title('Lyapunov exponents', fontsize=8)
        if has_resid and has_mode:
            # Stacked above (b) — hide x-label to avoid overlap
            ax_otd.set_xlabel('')
            ax_otd.tick_params(labelbottom=False)
        nk.panel_label(ax_otd, rf'$\bf{{({next(labels)})}}$')

    if has_resid and ax_res is not None:
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
        nk.panel_label(ax_res, rf'$\bf{{({next(labels)})}}$')

    for ff, ax_m, title in [
        (ff1, ax_mode1 if has_mode else None, 'OTD mode 1 ($v_y$)'),
        (ff2, ax_mode2 if has_mode and ff2 else None, 'OTD mode 2 ($v_y$)'),
    ]:
        if ff and ax_m is not None:
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
            nk.panel_label(ax_m, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
