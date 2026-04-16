#!/usr/bin/env python3
"""plot.py: Visualize Newton (dynamic + temperature) base flow results.

OUTPUTS: plot.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

CASE_DIR = Path(__file__).resolve().parent
OUTPUT = CASE_DIR / 'plot.png'


def make_white_anchored_cmap(name, colors):
    """Return a colormap whose zero level is rendered as white.

    The first color is forced to pure white so that quiescent regions are easy
    to separate from weak but nonzero structures in the wake.
    """
    return LinearSegmentedColormap.from_list(name, ['#ffffff', *colors])


VELOCITY_CMAP = LinearSegmentedColormap.from_list(
    'velocity_diverging_u1_white',
    ['#0b3c6d', '#4f8fc4', '#ffffff', '#f2b36f', '#b22222']
)

TEMPERATURE_CMAP = make_white_anchored_cmap(
    'temperature_white_heat',
    ['#fff3b0', '#fdc86d', '#f58b4c', '#d84b3a', '#7f0000']
)


def add_external_colorbar(fig, ax, mappable, ticks, tick_labels):
    """Place a horizontal colorbar above the axes, outside the plot area."""
    cbar = fig.colorbar(
        mappable,
        ax=ax,
        orientation='horizontal',
        location='top',
        fraction=0.06,
        pad=0.04
    )
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(tick_labels)
    cbar.outline.set_linewidth(0.5)
    return cbar


def main():
    nk.configure_style()

    has_conv = ((CASE_DIR / 'residu_newton.dat').exists() or
                (CASE_DIR / 'residu_arnoldi.dat').exists())
    bf_files = nk.find_fields('BF_*0.f*', CASE_DIR)
    if not bf_files:
        bf_files = nk.find_fields('BF*0.f*', CASE_DIR)

    ncols = int(has_conv) + 2  # velocity + temperature panels
    if not bf_files:
        ncols = max(int(has_conv), 1)

    ratios = []
    if has_conv:
        ratios.append(0.8)
    if bf_files:
        ratios.extend([1, 1])

    if not ratios:
        print('No data found — nothing to plot.')
        return

    fig = plt.figure(figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.52))
    gs = fig.add_gridspec(1, len(ratios), width_ratios=ratios, wspace=0.35)
    labels = iter('abcdefgh')
    col = 0

    if has_conv:
        ax_conv = fig.add_subplot(gs[col]); col += 1
        nk.plot_newton_convergence(ax_conv, CASE_DIR)
        handles, labels_conv = ax_conv.get_legend_handles_labels()
        if handles:
            ax_conv.legend(handles, labels_conv, fontsize=5, ncol=1,
                           handlelength=1.0, handletextpad=0.4,
                           borderpad=0.25, labelspacing=0.25,
                           columnspacing=0.6, loc='best')
        ax_conv.set_title('Convergence', fontsize=8)
        nk.panel_label(ax_conv, rf'$\bf{{({next(labels)})}}$')

    if bf_files:
        # Use the latest available saved base-flow-like snapshot in the folder.
        # This avoids plotting an older seed when both `BFre...` and `BF...`
        # files are present.
        x, y, fields, time = nk.read_field(bf_files[-1])
        triang = nk.make_triangulation(x, y)

        # Velocity magnitude
        ax_vel = fig.add_subplot(gs[col]); col += 1
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        vel_levels = np.linspace(0.0, 1.5, 257)
        vel_norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=1.5)
        cf = ax_vel.tricontourf(triang, umag, levels=vel_levels,
                                cmap=VELOCITY_CMAP, norm=vel_norm,
                                extend='both')
        add_external_colorbar(fig, ax_vel, cf,
                              ticks=[0.0, 1.0, 1.5],
                              tick_labels=['0', '1', '1.5'])
        ax_vel.set_aspect('equal')
        nk.add_cylinder_patches(ax_vel)
        ax_vel.set_xlim(-2, 20)
        ax_vel.set_ylim(-4, 4)
        ax_vel.set_xlabel(r'$x$', labelpad=-1)
        ax_vel.set_ylabel(r'$y$', labelpad=1)
        ax_vel.spines['right'].set_visible(False)
        ax_vel.spines['top'].set_visible(False)
        nk.panel_label(ax_vel, rf'$\bf{{({next(labels)})}}$')

        # Temperature
        ax_t = fig.add_subplot(gs[col])
        if 't' in fields:
            temp_max = max(1.0, float(np.nanmax(fields['t'])))
            cf2 = nk.tricontourf(ax_t, triang, fields['t'], levels=257,
                                 cmap=TEMPERATURE_CMAP, vmin=0.0,
                                 vmax=temp_max, extend='max')
            add_external_colorbar(fig, ax_t, cf2,
                                  ticks=[0.0, temp_max],
                                  tick_labels=['0', f'{temp_max:.1f}'])
        ax_t.set_aspect('equal')
        nk.add_cylinder_patches(ax_t)
        ax_t.set_xlim(-2, 20)
        ax_t.set_ylim(-4, 4)
        ax_t.set_xlabel(r'$x$', labelpad=-1)
        ax_t.spines['right'].set_visible(False)
        ax_t.spines['top'].set_visible(False)
        nk.panel_label(ax_t, rf'$\bf{{({next(labels)})}}$')

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close()


if __name__ == '__main__':
    main()
