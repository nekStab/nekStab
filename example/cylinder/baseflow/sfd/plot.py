#!/usr/bin/env python3
"""plot.py: Visualize SFD base flow results.

Convergence panel overlays SFD vs SFD_dyn to show that dynamic
parameters reach the same residual in fewer iterations (larger dt).

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
SFD_DYN_DIR = CASE_DIR.parent / 'sfd_dyn'
OUTPUT = CASE_DIR / 'plot.png'


def main():
    nk.configure_style()

    residu = CASE_DIR / 'residu.dat'
    residu_dyn = SFD_DYN_DIR / 'residu.dat'
    has_conv = residu.exists()
    has_dyn = residu_dyn.exists()
    bf_files = nk.find_fields('1cyl0.f*', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF_1cyl0.f*', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF*1cyl0.f*', CASE_DIR)

    if (has_conv or has_dyn) and bf_files:
        fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.45))
        gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.35)
        ax_conv = fig.add_subplot(gs[0])
        ax_bf = fig.add_subplot(gs[1])
    elif bf_files:
        fig, ax_bf = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.45))
        ax_conv = None
    else:
        fig, ax_conv = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5, nk.COL_WIDTH * 0.67))
        ax_bf = None

    labels = iter('abcdefgh')

    # ── Convergence: SFD vs SFD_dyn ──────────────────────────────────
    if ax_conv is not None and (has_conv or has_dyn):
        # Plot both residual histories against time (bottom axis)
        if has_conv:
            data_sfd = np.genfromtxt(str(residu))
            if data_sfd.ndim == 1:
                data_sfd = data_sfd.reshape(1, -1)
            t_sfd = data_sfd[:, 0]
            r_sfd = data_sfd[:, 1]
            ax_conv.semilogy(t_sfd, r_sfd, color='C0', lw=0.8,
                             label=f'SFD ({len(t_sfd)} iter.)')

        if has_dyn:
            data_dyn = np.genfromtxt(str(residu_dyn))
            if data_dyn.ndim == 1:
                data_dyn = data_dyn.reshape(1, -1)
            t_dyn = data_dyn[:, 0]
            r_dyn = data_dyn[:, 1]
            ax_conv.semilogy(t_dyn, r_dyn, color='C1', lw=0.8,
                             label=f'SFD dyn. ({len(t_dyn)} iter.)')

        ax_conv.set_xlabel(r'$t$')
        ax_conv.set_ylabel(r'$\|r\|$')
        ax_conv.legend(fontsize=5.5, loc='upper right')

        # ── Top axis: iteration count ────────────────────────────────
        ax_iter = ax_conv.twiny()
        t_min, t_max = ax_conv.get_xlim()

        # Tick marks based on SFD iteration spacing (finer grid = more iters)
        if has_conv and len(t_sfd) > 1:
            dt_sfd = (t_sfd[-1] - t_sfd[0]) / (len(t_sfd) - 1)
        else:
            dt_sfd = 1.0
        if has_dyn and len(t_dyn) > 1:
            dt_dyn = (t_dyn[-1] - t_dyn[0]) / (len(t_dyn) - 1)
        else:
            dt_dyn = dt_sfd

        # Show iteration ticks for the SFD case (finer steps)
        ax_iter.set_xlim(t_min / dt_sfd, t_max / dt_sfd)
        ax_iter.set_xlabel('iterations (SFD)', fontsize=6, labelpad=2)
        ax_iter.tick_params(labelsize=5)

        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

    # ── Base flow field ──────────────────────────────────────────────
    if ax_bf is not None and bf_files:
        x, y, fields, time = nk.read_field(bf_files[0])
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        cf = nk.tricontourf(ax_bf, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                          width="50%", height="5%", loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax_bf.set_aspect('equal')
        nk.add_cylinder_patches(ax_bf)
        ax_bf.set_xlim(-2, 20)
        ax_bf.set_ylim(-4, 4)
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
