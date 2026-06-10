#!/usr/bin/env python3
"""plot.py: Cylinder Re=180 adjoint Floquet stability — canonical stability panels.

Panel sequence (stability lane: BF -> SPECTRUM -> MODE):
  plot_bf.png       - UPO base flow (BF_*0.f00001 if present, else skipped)
  plot_spectrum.png - Floquet spectrum on unit circle (Spectre_Ha.dat)
  plot_mode.png     - adjoint Floquet mode (aRe1cyl0.f00001, vy field, RdBu)

Geometry: cylinder, xlim=(-2,20), ylim=(-4,4).

OUTPUTS: plot_spectrum.png, plot_mode.png [, plot_bf.png]
USAGE:   python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent


def find_bf(case_dir):
    for pat in ('BF_*0.f00001', 'BF*0.f00001'):
        matches = sorted(case_dir.glob(pat))
        if matches:
            return matches[-1]
    return None


def find_mode(case_dir):
    """Find adjoint Floquet mode real part."""
    matches = sorted(case_dir.glob('aRe*0.f00001'))
    if not matches:
        matches = sorted(case_dir.glob('aRe*0.f0*'))
    return matches[0] if matches else None


def render_cyl_field(ax, field_path, what='umag', cmap='Blues',
                     vmin=0, vmax=1.5, extend='max', symmetric=False):
    x, y, fields, _ = nk.read_field(str(field_path))
    triang = nk.make_triangulation(x, y)
    if what == 'umag':
        q = np.sqrt(fields['vx']**2 + fields['vy']**2)
    else:
        q = fields.get(what, fields.get('vy', fields.get('vx')))
    if q is None:
        return None
    if symmetric:
        bd = float(np.nanpercentile(np.abs(q), 99))
        vmin, vmax = -bd, bd
        extend = 'both'
    cf = nk.tricontourf(ax, triang, q, levels=257,
                        cmap=cmap, vmin=vmin, vmax=vmax, extend=extend)
    nk.add_cylinder_patches(ax)
    ax.set_xlim(-2, 20)
    ax.set_ylim(-4, 4)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$', labelpad=-1)
    ax.set_ylabel(r'$y$', labelpad=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return cf


def main():
    nk.configure_style()

    # PANEL 1 — BF (UPO base flow, skipped gracefully if not present)
    bf = find_bf(CASE_DIR)
    if bf:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        cf = render_cyl_field(ax, bf, what='umag', cmap='Blues',
                              vmin=0, vmax=1.5, extend='max')
        if cf is not None:
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width='50%', height='5%', loc=9,
                              ticks=[0, 1.5], tick_labels=['0', '1.5'])
        fig.savefig(CASE_DIR / 'plot_bf.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_bf.png"}')
        plt.close(fig)
    else:
        print('INFO: no BF field found — plot_bf.png skipped')

    # PANEL 2 — SPECTRUM (adjoint Floquet unit circle)
    spec_path = CASE_DIR / 'Spectre_Ha.dat'
    if not spec_path.exists():
        spec_path = CASE_DIR / 'Spectre_NSa.dat'
    if spec_path.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.8, nk.COL_WIDTH * 0.8))
        nk.setup_unit_circle_axes(ax, lim=1.5)
        nk.plot_spectrum_H_paper(ax, spec_path, label=r'$Re=180$')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.savefig(CASE_DIR / 'plot_spectrum.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_spectrum.png"}')
        plt.close(fig)

    # PANEL 3 — MODE (adjoint Floquet vy, RdBu)
    mode = find_mode(CASE_DIR)
    if mode:
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
        cf = render_cyl_field(ax, mode, what='vy', cmap='RdBu',
                              symmetric=True)
        if cf is not None:
            bd = float(max(abs(cf.get_clim()[0]), abs(cf.get_clim()[1])))
            bdr = round(bd, 2)
            nk.inset_colorbar(ax, cf, orientation='horizontal',
                              width='50%', height='5%', loc=9,
                              ticks=[-bdr, bdr],
                              tick_labels=[f'{-bdr}', f'{bdr}'])
        fig.savefig(CASE_DIR / 'plot_mode.png', dpi=600, bbox_inches='tight')
        print(f'Saved {CASE_DIR / "plot_mode.png"}')
        plt.close(fig)
    else:
        print('INFO: no adjoint Floquet mode file found — plot_mode.png skipped')


if __name__ == '__main__':
    main()
