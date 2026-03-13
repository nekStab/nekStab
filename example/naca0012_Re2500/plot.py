#!/usr/bin/env python3
"""plot.py: Visualize NACA 0012 Re=2500 stability results.

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

    spec_ns = CASE_DIR / 'Spectre_NSd.dat'
    bf_files = nk.find_fields('BF_naca00120.f*', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('dns_*naca00120.f*', CASE_DIR)

    has_spec = spec_ns.exists()
    has_bf = bool(bf_files)

    ncols = int(has_spec) + int(has_bf)
    ncols = max(ncols, 1)

    ratios = []
    if has_spec:
        ratios.append(0.8)
    if has_bf:
        ratios.append(1)
    if not ratios:
        ratios = [1]

    fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.5))
    gs = fig.add_gridspec(1, len(ratios), width_ratios=ratios, wspace=0.35)
    labels = iter('abcdefgh')
    col = 0

    if has_spec:
        ax_spec = fig.add_subplot(gs[col]); col += 1
        nk.plot_spectrum_NS_paper(ax_spec, spec_ns, label=r'$Re=2500$')
        ylim = ax_spec.get_ylim()
        ax_spec.axhspan(min(ylim[0], -0.2), 0, facecolor='gray', alpha=0.3, zorder=-1)
        ax_spec.set_title(r'$\sigma$ vs $f$', fontsize=8)
        nk.panel_label(ax_spec, rf'$\bf{{({next(labels)})}}$')

    if has_bf:
        ax_bf = fig.add_subplot(gs[col])
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        cf = nk.tricontourf(ax_bf, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax_bf.set_aspect('equal')
        nk.add_naca0012_patch(ax_bf)
        ax_bf.set_xlim(-1, 5)
        ax_bf.set_ylim(-2, 2)
        ax_bf.set_xlabel(r'$x$', labelpad=-1)
        ax_bf.set_ylabel(r'$y$', labelpad=1)
        ax_bf.spines['right'].set_visible(False)
        ax_bf.spines['top'].set_visible(False)
        nk.panel_label(ax_bf, rf'$\bf{{({next(labels)})}}$')

    if not has_spec and not has_bf:
        ax = fig.add_subplot(gs[0])
        ax.text(0.5, 0.5, 'No data found', transform=ax.transAxes, ha='center')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
