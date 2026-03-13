#!/usr/bin/env python3
"""plot.py: Visualize direct stability analysis results.

OUTPUTS: plot.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
OUTPUT = CASE_DIR / 'plot.png'


def main():
    nk.configure_style()

    spec_ns = CASE_DIR / 'Spectre_NSd.dat'
    dRe_files = sorted(CASE_DIR.glob('dRe*0.f0*'))

    has_spec = spec_ns.exists()
    has_mode = bool(dRe_files)

    if has_spec and has_mode:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_spec = fig.add_subplot(gs[0])
        ax_mode = fig.add_subplot(gs[1])
    elif has_mode:
        fig, ax_mode = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.45))
        ax_spec = None
    else:
        fig, ax_spec = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5, nk.COL_WIDTH * 0.67))
        ax_mode = None

    labels = iter('abcdefgh')

    if ax_spec is not None and has_spec:
        nk.plot_spectrum_NS_paper(ax_spec, spec_ns, label=r'$Re=50$')
        ylim = ax_spec.get_ylim()
        ax_spec.axhspan(min(ylim[0], -0.2), 0, facecolor='gray', alpha=0.3, zorder=-1)
        ax_spec.set_title(r'$\sigma$ vs $f$', fontsize=8)
        nk.panel_label(ax_spec, rf'$\bf{{({next(labels)})}}$')

    if ax_mode is not None and has_mode:
        x, y, fields, time = nk.read_field(dRe_files[0])
        triang = nk.make_triangulation(x, y)
        q = fields.get('vy', fields.get('vx'))
        if q is not None:
            bd = np.nanpercentile(np.abs(q), 99)
            cf = nk.tricontourf(ax_mode, triang, q, levels=257,
                                cmap='RdBu', vmin=-bd, vmax=bd, extend='both')
            bdr = round(bd, 2)
            nk.inset_colorbar(ax_mode, cf, orientation='horizontal',
                              width="50%", height="5%", loc=9,
                              ticks=[-bdr, bdr],
                              tick_labels=[f'{-bdr}', f'{bdr}'])
        ax_mode.set_aspect('equal')
        nk.add_cylinder_patches(ax_mode)
        ax_mode.set_xlim(-2, 20)
        ax_mode.set_ylim(-4, 4)
        ax_mode.set_xlabel(r'$x$', labelpad=-1)
        ax_mode.set_ylabel(r'$y$', labelpad=1)
        ax_mode.spines['right'].set_visible(False)
        ax_mode.spines['top'].set_visible(False)
        nk.panel_label(ax_mode, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
