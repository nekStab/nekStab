#!/usr/bin/env python3
"""plot.py: NACA 0012 Re=2000 Newton base flow — canonical baseflow panels.

Panel sequence (canonical baseflow lane: IC -> RESIDUAL -> BF):
  plot_ic.png       - initial condition (seed file named by startFrom, |u|, Blues)
  plot_residual.png - Newton convergence history
  plot_bf.png       - converged base flow (|u|, Blues)

Geometry: NACA 0012 airfoil, xlim=(-1,5), ylim=(-2,2).

OUTPUTS: plot_ic.png, plot_residual.png, plot_bf.png
USAGE:   python plot.py
"""
from pathlib import Path
import re
import sys
sys.path.insert(0, str(next(d for d in Path(__file__).resolve().parents if (d / 'nekplot.py').exists())))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent

XLIM = (-1, 5)
YLIM = (-2, 2)


def find_snapshot(case_dir):
    """Return the base flow that Newton wrote. The seed BF_seed_* does not match these patterns."""
    for pat in ('BF_naca00120.f00001', 'BF_naca00120.f*'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            return matches[-1]
    return None


def find_seed(case_dir):
    """Return the file named by startFrom in the .par file, or None."""
    for par in sorted(case_dir.glob('*.par')):
        found = re.search(r'^\s*startFrom\s*=\s*([^\s#]+)', par.read_text(), re.M)
        if found and (case_dir / found.group(1)).exists():
            return case_dir / found.group(1)
    return None


def render_naca_field(ax, field_path):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
    cf = nk.tricontourf(ax, triang, umag, levels=257,
                        cmap='Blues', vmin=0, vmax=1.5, extend='max')
    nk.inset_colorbar(ax, cf, orientation='horizontal',
                      width='50%', height='5%', loc=9,
                      ticks=[0, 1.5], tick_labels=['0', '1.5'])
    nk.add_naca0012_patch(ax)
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def main():
    nk.configure_style()

    # PANEL 1 — IC (the DNS snapshot that Newton starts from)
    seed = find_seed(CASE_DIR)
    if seed:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.5))
        render_naca_field(ax, seed)
        fig.savefig(CASE_DIR / 'plot_ic.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_ic.png"}')
        plt.close(fig)

    # PANEL 2 — RESIDUAL (Newton convergence)
    if (CASE_DIR / 'residu_newton.dat').exists() or (CASE_DIR / 'residu_arnoldi.dat').exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        nk.plot_newton_convergence(ax, CASE_DIR)
        handles, labs = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labs, fontsize=6, loc='best')
        ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_residual.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_residual.png"}')
        plt.close(fig)

    # PANEL 3 — BF (converged base flow)
    bf = find_snapshot(CASE_DIR)
    if bf:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.5))
        render_naca_field(ax, bf)
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)
    else:
        print('WARNING: no base-flow field found — plot_bf.png skipped')


if __name__ == '__main__':
    main()
