#!/usr/bin/env python3
"""plot.py: Cylinder Re=30 thermal Newton (dynamic + temperature) base flow.

Panel sequence:
  plot_ic.png           - initial condition (seed from startFrom, if exists)
  plot_residual.png     - Newton convergence history (renamed from plot_convergence)
  plot_bf_velocity.png  - converged BF velocity magnitude (if BF field exists)
  plot_bf_temperature.png - converged BF temperature field (if BF field has 't')

Note: BF field for the current 2128-element mesh is not yet available (seed
archived at .archive_mesh_1996_2026-05-21/). plot_bf_*.png panels are skipped
gracefully when no BF file is present.

Geometry: cylinder, xlim=(-2,20), ylim=(-4,4).

OUTPUTS: [plot_ic.png,] plot_residual.png [, plot_bf_velocity.png, plot_bf_temperature.png]
USAGE:   python plot.py
"""
from pathlib import Path
import sys, re
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

CASE_DIR = Path(__file__).resolve().parent

VELOCITY_CMAP = LinearSegmentedColormap.from_list(
    'velocity_diverging_u1_white',
    ['#0b3c6d', '#4f8fc4', '#ffffff', '#f2b36f', '#b22222']
)

TEMPERATURE_CMAP = LinearSegmentedColormap.from_list(
    'temperature_white_heat',
    ['#ffffff', '#fff3b0', '#fdc86d', '#f58b4c', '#d84b3a', '#7f0000']
)


def parse_startfrom(par_path):
    if not par_path.exists():
        return None
    txt = par_path.read_text()
    m = re.search(r'^\s*startFrom\s*=\s*([^\s#]+)', txt, re.M | re.I)
    if not m:
        return None
    val = m.group(1).strip().strip('"').strip("'")
    if val in ('0', ''):
        return None
    return val


def find_ic(case_dir):
    """Find IC seed field. Skips commented-out entries in par file."""
    par_files = list(case_dir.glob('*.par'))
    if not par_files:
        return None
    ic_name = parse_startfrom(par_files[0])
    if not ic_name:
        return None
    for d in (case_dir, case_dir.parent, case_dir.parent.parent):
        p = d / ic_name
        if p.exists():
            return p
    return None


def find_bf(case_dir):
    """Find converged BF field (excludes archived incompatible meshes)."""
    for pat in ('BF_*0.f00001', 'BF*0.f00001'):
        matches = [f for f in sorted(case_dir.glob(pat))
                   if '.archive' not in str(f)]
        if matches:
            return matches[-1]
    return None


def add_external_colorbar(fig, ax, mappable, ticks, tick_labels):
    """Place a horizontal colorbar above the axes, outside the plot area."""
    cbar = fig.colorbar(
        mappable, ax=ax, orientation='horizontal', location='top',
        fraction=0.06, pad=0.04
    )
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(tick_labels)
    cbar.outline.set_linewidth(0.5)
    return cbar


def main():
    nk.configure_style()

    # PANEL 1 — IC (seed from startFrom; skipped gracefully if none exists)
    ic = find_ic(CASE_DIR)
    if ic:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        x, y, fields, _ = nk.read_field(str(ic))
        triang = nk.make_triangulation(x, y)
        if 'vx' in fields and 'vy' in fields:
            umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
            vel_levels = np.linspace(0.0, 1.5, 257)
            cf = ax.tricontourf(triang, umag, levels=vel_levels,
                                cmap='Blues', extend='max')
            add_external_colorbar(fig, ax, cf,
                                  ticks=[0.0, 1.5],
                                  tick_labels=['0', '1.5'])
        nk.add_cylinder_patches(ax)
        ax.set_xlim(-2, 20)
        ax.set_ylim(-4, 4)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$x$', labelpad=-1)
        ax.set_ylabel(r'$y$', labelpad=1)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_ic.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_ic.png"}')
        plt.close(fig)
    else:
        print('INFO: no IC seed found (startFrom commented out) — plot_ic.png skipped')

    # PANEL 2 — RESIDUAL (Newton convergence)
    residu_nwt = CASE_DIR / 'residu_newton.dat'
    residu_arn = CASE_DIR / 'residu_arnoldi.dat'
    if residu_nwt.exists() or residu_arn.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_newton_convergence(ax, CASE_DIR)
        handles, labs_conv = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labs_conv, fontsize=5, ncol=1,
                      handlelength=1.0, handletextpad=0.4,
                      borderpad=0.25, labelspacing=0.25,
                      columnspacing=0.6, loc='best')
        ax.set_title('Convergence', fontsize=8)
        ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_residual.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_residual.png"}')
        plt.close(fig)

    # PANELS 3 & 4 — BF velocity and temperature (skipped if no BF file)
    bf = find_bf(CASE_DIR)
    if bf:
        x, y, fields, time = nk.read_field(str(bf))
        triang = nk.make_triangulation(x, y)

        # BF velocity magnitude
        if 'vx' in fields and 'vy' in fields:
            fig, ax_vel = plt.subplots(1, 1,
                                       figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
            umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
            vel_levels = np.linspace(0.0, 1.5, 257)
            vel_norm = TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=1.5)
            cf = ax_vel.tricontourf(triang, umag, levels=vel_levels,
                                    cmap=VELOCITY_CMAP, norm=vel_norm,
                                    extend='both')
            add_external_colorbar(fig, ax_vel, cf,
                                  ticks=[0.0, 1.0, 1.5],
                                  tick_labels=['0', '1', '1.5'])
            nk.add_cylinder_patches(ax_vel)
            ax_vel.set_xlim(-2, 20)
            ax_vel.set_ylim(-4, 4)
            ax_vel.set_aspect('equal')
            ax_vel.set_xlabel(r'$x$', labelpad=-1)
            ax_vel.set_ylabel(r'$y$', labelpad=1)
            ax_vel.spines['right'].set_visible(False)
            ax_vel.spines['top'].set_visible(False)
            fig.savefig(CASE_DIR / 'plot_bf_velocity.png', dpi=600, bbox_inches='tight')
            print(f'Saved {CASE_DIR / "plot_bf_velocity.png"}')
            plt.close(fig)

        # BF temperature
        if 't' in fields:
            fig, ax_t = plt.subplots(1, 1,
                                     figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
            temp_max = max(1.0, float(np.nanmax(fields['t'])))
            cf2 = nk.tricontourf(ax_t, triang, fields['t'], levels=257,
                                 cmap=TEMPERATURE_CMAP, vmin=0.0,
                                 vmax=temp_max, extend='max')
            add_external_colorbar(fig, ax_t, cf2,
                                  ticks=[0.0, temp_max],
                                  tick_labels=['0', f'{temp_max:.1f}'])
            nk.add_cylinder_patches(ax_t)
            ax_t.set_xlim(-2, 20)
            ax_t.set_ylim(-4, 4)
            ax_t.set_aspect('equal')
            ax_t.set_xlabel(r'$x$', labelpad=-1)
            ax_t.set_ylabel(r'$y$', labelpad=1)
            ax_t.spines['right'].set_visible(False)
            ax_t.spines['top'].set_visible(False)
            fig.savefig(CASE_DIR / 'plot_bf_temperature.png', dpi=600, bbox_inches='tight')
            print(f'Saved {CASE_DIR / "plot_bf_temperature.png"}')
            plt.close(fig)
    else:
        print('INFO: no BF field file found — plot_bf_velocity.png and plot_bf_temperature.png skipped')


if __name__ == '__main__':
    main()
