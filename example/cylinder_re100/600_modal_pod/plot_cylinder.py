#!/usr/bin/env python3
"""
plot_cylinder.py: Visualize Nek5000 field files (velocity magnitude).

OUTPUTS: plot_cylinder.png
USAGE:   p plot_cylinder.py [options] [file]

Options:
    --stream     Overlay streamlines
    --xlim X0,X1 X-axis limits (default: -2,20)
    --ylim Y0,Y1 Y-axis limits (default: -4,4)

Examples:
    p plot_cylinder.py                          # velocity, first snapshot
    p plot_cylinder.py --stream 1cyl0.f00050    # with streamlines
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from pymech.neksuite import readnek
from scipy.interpolate import griddata

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DPI = 300
FORMAT = 'png'
CMAP_VEL = 'Blues'

SCRIPT = Path(__file__).resolve()
OUTPUT_PNG = SCRIPT.with_suffix('.png')


# ---------------------------------------------------------------------------
# Data extraction
# ---------------------------------------------------------------------------

def extract_raw(field):
    """Flatten element arrays into 1-D coordinate/velocity vectors."""
    nel = field.nel
    nx, ny, nz = field.lr1
    npe = nx * ny * max(nz, 1)
    ntotal = nel * npe

    x  = np.empty(ntotal)
    y  = np.empty(ntotal)
    ux = np.empty(ntotal)
    uy = np.empty(ntotal)

    for ie, elem in enumerate(field.elem):
        s = slice(ie * npe, (ie + 1) * npe)
        x[s]  = elem.pos[0].ravel()
        y[s]  = elem.pos[1].ravel()
        ux[s] = elem.vel[0].ravel()
        uy[s] = elem.vel[1].ravel()

    return x, y, ux, uy


def interpolate_to_grid(x, y, fields, xlim, ylim, ngrid=400):
    """Interpolate scattered SEM data to a uniform grid.

    Uses cubic interpolation with linear/nearest fallback for points
    near domain boundaries (e.g. cylinder wall) where cubic fails.

    Parameters
    ----------
    x, y   : 1-D arrays of SEM nodal coordinates
    fields : dict of {name: values} to interpolate
    xlim, ylim : plot limits (determines grid extent)
    ngrid  : number of grid points along x

    Returns
    -------
    xi, yi : 1-D grid vectors
    result : dict of {name: 2-D interpolated array}
    """
    aspect = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])
    nxi = ngrid
    nyi = max(int(ngrid * aspect), 50)

    xi = np.linspace(xlim[0], xlim[1], nxi)
    yi = np.linspace(ylim[0], ylim[1], nyi)
    Xi, Yi = np.meshgrid(xi, yi)

    pts = np.column_stack([x, y])
    result = {}
    for name, vals in fields.items():
        Z = griddata(pts, vals, (Xi, Yi), method='cubic')
        # Fill near-boundary NaN with linear, then nearest
        mask = np.isnan(Z)
        if mask.any():
            Z_lin = griddata(pts, vals, (Xi, Yi), method='linear')
            Z[mask] = Z_lin[mask]
            mask = np.isnan(Z)
            if mask.any():
                Z_near = griddata(pts, vals, (Xi, Yi), method='nearest')
                Z[mask] = Z_near[mask]
        result[name] = Z

    return xi, yi, Xi, Yi, result


def mask_body(Xi, Yi, *arrays, radius=0.51):
    """Mask grid points inside the cylinder body."""
    r = np.sqrt(Xi**2 + Yi**2)
    return [np.ma.masked_where(r < radius, a) for a in arrays]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_field(field, output_path=None, show_stream=False,
               xlim=(-2, 20), ylim=(-4, 4)):
    """Render velocity magnitude from a Nek5000 snapshot."""
    if output_path is None:
        output_path = OUTPUT_PNG

    print(f'  Time: {field.time:.3f}')
    print(f'  Elements: {field.nel}, GLL: {field.lr1}')

    # --- extract ---
    print('  Extracting mesh data...')
    x, y, ux, uy = extract_raw(field)

    # --- interpolate to uniform grid ---
    fields = {'ux': ux, 'uy': uy}

    print('  Interpolating to uniform grid...')
    xi, yi, Xi, Yi, interp = interpolate_to_grid(x, y, fields, xlim, ylim)

    Ux = interp['ux']
    Uy = interp['uy']
    vel_mag = np.sqrt(Ux**2 + Uy**2)

    # Mask body interior
    [vel_mag_m] = mask_body(Xi, Yi, vel_mag)

    # --- helpers ---
    theta = np.linspace(0, 2 * np.pi, 100)
    cyl_x, cyl_y = 0.5 * np.cos(theta), 0.5 * np.sin(theta)

    def add_cylinder(ax):
        ax.fill(cyl_x, cyl_y, 'white', ec='k', lw=1.0, zorder=10)

    def add_streamlines(ax):
        Ux_s, Uy_s = Ux.copy(), Uy.copy()
        r = np.sqrt(Xi**2 + Yi**2)
        Ux_s = np.ma.masked_where(r < 0.55, Ux_s)
        Uy_s = np.ma.masked_where(r < 0.55, Uy_s)
        ax.streamplot(xi, yi, Ux_s, Uy_s, color='k', linewidth=0.4,
                      density=1.2, arrowsize=0.5, arrowstyle='->')

    # --- figure ---
    fig, ax = plt.subplots(figsize=(12, 4))
    vmax_vel = min(1.5, np.nanpercentile(vel_mag, 99.5))
    tcf = ax.pcolormesh(xi, yi, vel_mag_m, shading='gouraud',
                        cmap=CMAP_VEL, vmin=0, vmax=vmax_vel)
    if show_stream:
        add_streamlines(ax)
    fig.colorbar(tcf, ax=ax, shrink=0.8, pad=0.02,
                 label=r'$|\mathbf{u}|$')
    add_cylinder(ax)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Velocity magnitude, t = {field.time:.2f}')

    fig.tight_layout()
    fig.savefig(output_path, dpi=DPI, bbox_inches='tight')
    print(f'  Saved {output_path}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    args = sys.argv[1:]
    show_stream = False
    xlim = (-2, 20)
    ylim = (-4, 4)
    fname = None

    i = 0
    while i < len(args):
        a = args[i]
        if   a == '--stream': show_stream = True
        elif a == '--xlim' and i + 1 < len(args):
            i += 1; x0, x1 = map(float, args[i].split(',')); xlim = (x0, x1)
        elif a == '--ylim' and i + 1 < len(args):
            i += 1; y0, y1 = map(float, args[i].split(',')); ylim = (y0, y1)
        elif not a.startswith('-') and not a.endswith('.py'):
            fname = a
        i += 1

    if fname is None:
        files = sorted(Path('.').glob('*0.f?????'))
        if not files:
            print('No field files found (*0.f?????)')
            sys.exit(1)
        fname = str(files[0])

    return fname, show_stream, xlim, ylim


def main():
    fname, show_stream, xlim, ylim = parse_args()
    print(f'Reading {fname}...')
    field = readnek(fname)
    plot_field(field, show_stream=show_stream, xlim=xlim, ylim=ylim)


if __name__ == '__main__':
    main()
