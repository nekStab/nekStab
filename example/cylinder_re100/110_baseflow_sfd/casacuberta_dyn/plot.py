#!/usr/bin/env python3
"""plot.py: Cylinder SFD Casacuberta with dynamic tolerance — canonical baseflow panels.

Panel sequence (baseflow lane + extra dyn_tol panel):
  plot_ic.png                - initial condition (seed from startFrom)
  plot_residual.png          - SFD dyn vs fixed-tol residual, with auto-cap ylim
  plot_dynamic_tolerance.png - tolerance schedule (current/requested/used/cap)
  plot_bf.png                - converged base flow

OUTPUTS: plot_ic.png, plot_residual.png, plot_dynamic_tolerance.png, plot_bf.png
USAGE:   python plot.py
"""
from pathlib import Path
import sys, re
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import matplotlib.pyplot as plt
import numpy as np
import nekplot as nk

CASE_DIR = Path(__file__).resolve().parent
BASELINE_DIR = CASE_DIR.parent / 'casacuberta'  # SFD fixed-tol sibling for comparison


def parse_startfrom(par_path):
    if not par_path.exists(): return None
    m = re.search(r'^\s*startFrom\s*=\s*([^\s#]+)', par_path.read_text(), re.M)
    if not m: return None
    val = m.group(1).strip().strip('"').strip("'")
    return None if val in ('0', '') else val


def find_ic(case_dir):
    par_files = list(case_dir.glob('*.par'))
    if not par_files: return None
    name = parse_startfrom(par_files[0])
    if not name: return None
    for d in (case_dir, case_dir.parent, case_dir.parent.parent):
        p = d / name
        if p.exists(): return p
    return None


def find_bf(case_dir):
    for pat in ('BF_*0.f00001', 'BF*0.f00001', '1cyl0.f*'):
        m = sorted(case_dir.glob(pat))
        if m: return m[-1]
    return None


def render_cyl(ax, field_path, vmax=1.5):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    q = np.sqrt(fields['vx']**2 + fields['vy']**2)
    cf = nk.tricontourf(ax, triang, q, levels=257,
                        cmap='Blues', vmin=0, vmax=vmax, extend='max')
    nk.inset_colorbar(ax, cf, orientation='horizontal', width='50%',
                      height='5%', loc=9, ticks=[0, vmax],
                      tick_labels=['0', f'{vmax}'])
    ax.set_aspect('equal')
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 20); ax.set_ylim(-4, 4)
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def main():
    nk.configure_style()

    # PANEL 1 — IC
    ic = find_ic(CASE_DIR)
    if ic:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        render_cyl(ax, ic)
        fig.savefig(CASE_DIR / 'plot_ic.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_ic.png"}')
        plt.close(fig)

    # PANEL 2 — RESIDUAL (SFD dyn vs fixed-tol baseline, auto-ylim from nk)
    residu = CASE_DIR / 'residu.dat'
    baseline_residu = BASELINE_DIR / 'residu.dat'
    if residu.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        if baseline_residu.exists():
            nk.plot_residuals(ax, baseline_residu, label='Casacuberta fixed', color='0.45')
        nk.plot_residuals(ax, residu, label='Casacuberta dyn', color='C1')
        ax.axhline(1e-9, color='k', ls=':', lw=0.5, label=r'tol $10^{-9}$')
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r'$\|r\|$')
        ax.legend(fontsize=7, loc='upper right')
        ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False); ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_residual.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_residual.png"}')
        plt.close(fig)

    # PANEL 3 — DYNAMIC TOLERANCE schedule (extra panel for this lane)
    dyn_tol = CASE_DIR / 'dyn_tol.dat'
    if residu.exists() and dyn_tol.exists():
        residu_data = np.genfromtxt(str(residu))
        dyn_data = np.genfromtxt(str(dyn_tol), comments='#')
        if residu_data.ndim == 1: residu_data = residu_data.reshape(1, -1)
        if dyn_data.ndim == 1: dyn_data = dyn_data.reshape(1, -1)

        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH, nk.COL_WIDTH * 0.67))
        if residu_data.shape[1] >= 4:
            ax.semilogy(residu_data[:, 0], residu_data[:, 3], color='0.35',
                        lw=0.7, label='current tol')
        if dyn_data.shape[1] >= 4:
            ax.semilogy(dyn_data[:, 0], dyn_data[:, 3], 'o-', ms=2.0,
                        lw=0.6, color='C1', label='requested')
        if dyn_data.shape[1] >= 5:
            ax.semilogy(dyn_data[:, 0], dyn_data[:, 4], 'o-', ms=2.0,
                        lw=0.6, color='C2', label='used')
        if dyn_data.shape[1] >= 6 and np.any(dyn_data[:, 5] > 0.0):
            ax.semilogy(dyn_data[:, 0], dyn_data[:, 5], '--',
                        lw=0.6, color='C3', label='cap')
        ax.set_xlabel(r'$t$')
        ax.set_ylabel('solver tol')
        ax.legend(fontsize=6, loc='best')
        ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)
        ax.spines['right'].set_visible(False); ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_dynamic_tolerance.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_dynamic_tolerance.png"}')
        plt.close(fig)

    # PANEL 4 — BF
    bf = find_bf(CASE_DIR)
    if bf:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        render_cyl(ax, bf)
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)


if __name__ == '__main__':
    main()
