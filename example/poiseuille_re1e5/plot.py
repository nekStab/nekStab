#!/usr/bin/env python3
"""plot.py: Visualize RANS Poiseuille channel stability results.

Panel (a): SFD convergence (residu.dat) or NS spectrum (Spectre_NSd.dat)
Panel (b): Base flow velocity magnitude (Blues colormap)

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

    residu = CASE_DIR / 'residu.dat'
    spec_ns = CASE_DIR / 'Spectre_NSd.dat'
    bf_files = nk.find_fields('BF_*0.f*', CASE_DIR)

    has_conv = residu.exists()
    has_spec = spec_ns.exists()
    has_bf = bool(bf_files)
    has_left = has_conv or has_spec

    if has_left and has_bf:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_left = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif has_bf:
        fig, ax_bf = plt.subplots(1, 1,
                                  figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.45))
        ax_left = None
    elif has_left:
        fig, ax_left = plt.subplots(1, 1,
                                    figsize=(nk.COL_WIDTH * 0.5,
                                             nk.COL_WIDTH * 0.67))
        ax_bf = None
    else:
        print('No data files found — nothing to plot.')
        return

    labels = iter('abcdefgh')

    # ── Panel (a): convergence or spectrum ──
    if ax_left is not None:
        if has_conv:
            nk.plot_residuals(ax_left, residu)
            ax_left.set_title('SFD convergence', fontsize=8)
        elif has_spec:
            nk.plot_spectrum_NS_paper(ax_left, spec_ns,
                                      label=r'$Re=10^5$')
            ylim = ax_left.get_ylim()
            ax_left.axhspan(min(ylim[0], -0.2), 0,
                            facecolor='gray', alpha=0.3, zorder=-1)
            ax_left.set_title(r'$\sigma$ vs $f$', fontsize=8)
        nk.panel_label(ax_left, rf'$\bf{{({next(labels)})}}$')

    # ── Panel (b): base flow velocity magnitude ──
    if ax_bf is not None and has_bf:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        bd = np.nanpercentile(umag, 99.5)
        cf = nk.tricontourf(ax_bf, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=bd, extend='max')
        bdr = round(bd, 2)
        nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, bdr],
                          tick_labels=['0', f'{bdr}'])
        ax_bf.set_aspect('equal')
        ax_bf.set_xlabel(r'$x$', labelpad=-1)
        ax_bf.set_ylabel(r'$y$', labelpad=1)
        ax_bf.spines['right'].set_visible(False)
        ax_bf.spines['top'].set_visible(False)
        nk.panel_label(ax_bf, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
