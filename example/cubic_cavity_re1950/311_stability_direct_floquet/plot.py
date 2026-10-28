#!/usr/bin/env python3
"""plot.py: Visualize cubic-cavity direct Floquet results (Re=1950, 3D).

Floquet multiplier spectrum on the unit circle (geometry-independent).
The leading mode is 3D and is not drawn: the shared nekplot helpers make 2D contour plots, and a 2D
projection of a 3D mode would mislead. The mode is in dRecav0.f00001 and dImcav0.f00001.

The base flow is the on-orbit snapshot BF_cav0.f00001 (t = 1700): the speed |u| on the planes z = 0 (a) and
x = 0 (b), as the other stages of this case draw it (the cavity spans [-0.5, 0.5]^3 in units of the edge L; the
planes are faces of the elements, so their nodes are exact).

OUTPUTS: plot_spectrum.png, plot_bf.png
USAGE:   PYTHONPATH=<repo>/example python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt
import numpy as np
from pymech.neksuite import readnek

CASE_DIR = Path(__file__).resolve().parent


def plane_speed(path, normal):
    """Return the two in-plane coordinates and |u| on the plane (normal = 0 or 2: x = 0 or z = 0).

    A node on an element edge belongs to several elements; its entries are averaged.
    """
    data = readnek(str(path))
    pos = np.concatenate([e.pos.reshape(3, -1) for e in data.elem], axis=1)
    speed = np.concatenate([np.sqrt(np.sum(e.vel**2, axis=0)).ravel() for e in data.elem])
    on_plane = np.isclose(pos[normal], 0.0, atol=1e-9)
    in_plane = [a for a in (0, 1, 2) if a != normal]
    keys = np.round(pos[in_plane][:, on_plane].T, 9)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = inverse.ravel()
    average = np.bincount(inverse, weights=speed[on_plane]) / np.bincount(inverse)
    return unique[:, 0], unique[:, 1], average


def plot_baseflow(path, out):
    """Draw |u| of the orbit snapshot on the planes z = 0 (a) and x = 0 (b)."""
    fig, axes = plt.subplots(1, 2, figsize=(nk.COL_WIDTH * 1.8, nk.COL_WIDTH * 1.0))
    for ax, normal, tag, names in ((axes[0], 2, 'a', ('x', 'y')), (axes[1], 0, 'b', ('y', 'z'))):
        u, v, speed = plane_speed(path, normal)
        top = round(float(np.percentile(speed, 99.5)), 2)
        cf = ax.tricontourf(u, v, speed, levels=np.linspace(0.0, top, 25), cmap='viridis', extend='max')
        ax.set_aspect('equal')
        ax.set_xlabel(rf'${names[0]}$ ($L$)')
        ax.set_ylabel(rf'${names[1]}$ ($L$)')
        ax.set_title(rf'({tag}) $|\mathbf{{u}}|$ at ${"zx"[tag == "b"]}=0$', loc='left', fontsize=7)
        cb = fig.colorbar(cf, ax=ax, orientation='horizontal', shrink=0.9, pad=0.12,
                          ticks=np.linspace(0.0, top, 3), format='%.2f')
        cb.set_label(r'$|\mathbf{u}|$ ($U$)')
    fig.savefig(out, dpi=400, bbox_inches='tight')
    print(f'Saved {out}')
    plt.close(fig)


def main():
    nk.configure_style()

    bf = CASE_DIR / 'BF_cav0.f00001'
    if bf.exists():
        plot_baseflow(bf, CASE_DIR / 'plot_bf.png')

    spec_h = CASE_DIR / 'Spectre_Hd.dat'
    if spec_h.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                              nk.COL_WIDTH * 0.5))
        nk.setup_unit_circle_axes(ax, lim=1.5)
        nk.plot_spectrum_H_paper(ax, spec_h, label=r'$Re = 1950$')
        out = CASE_DIR / 'plot_spectrum.png'
        fig.savefig(out, dpi=600, bbox_inches='tight')
        print(f'Saved {out}')
        plt.close(fig)
    else:
        print('Spectre_Hd.dat not found; nothing to plot')


if __name__ == '__main__':
    main()
