#!/usr/bin/env python3
"""
plot_modes.py: Visualize spatial modes from Nek5000 modal analysis.

OUTPUTS: modes_mean.png, modes_pod.png, modes_dmd.png, modes_spod.png
USAGE:   p plot_modes.py              # all decompositions, 2 modes each
         p plot_modes.py -n 3         # 3 modes each
         p plot_modes.py -n 2 --spod  # SPOD only, 2 modes

Options:
    -n N             Number of modes to plot per decomposition (default: 2)
    --pod/--dmd/--spod/--mean   Select specific decomposition(s)
    --xlim X0,X1     X-axis limits (default: auto from mesh)
    --ylim Y0,Y1     Y-axis limits (default: auto from mesh)
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
CMAP = 'RdBu_r'
CMAP_SEQ = 'Blues'
SCRIPT = Path(__file__).resolve()
HERE = SCRIPT.parent


# ---------------------------------------------------------------------------
# Core I/O
# ---------------------------------------------------------------------------

def extract_raw(field):
    """Flatten spectral element data into 1-D arrays."""
    nel = field.nel
    nx, ny, nz = field.lr1
    npts = nx * ny * max(nz, 1)
    ntotal = nel * npts

    x  = np.empty(ntotal)
    y  = np.empty(ntotal)
    ux = np.empty(ntotal)
    uy = np.empty(ntotal)

    for ie, elem in enumerate(field.elem):
        s = slice(ie * npts, (ie + 1) * npts)
        x[s]  = elem.pos[0].ravel()
        y[s]  = elem.pos[1].ravel()
        ux[s] = elem.vel[0].ravel()
        uy[s] = elem.vel[1].ravel()

    return x, y, ux, uy


def find_mode_files(prefix):
    """Find Nek5000 field files matching a prefix, sorted by index."""
    return sorted(HERE.glob(f'{prefix}*0.f?????'))


def auto_limits(x, y):
    """Compute view limits focused on the near-wake region.

    Detects the body (minimum-radius cluster) and sets limits to show
    the near-wake: ~4 diameters upstream, ~40 diameters downstream,
    ±6 diameters in the transverse direction.
    """
    r = np.sqrt(x**2 + y**2)
    r_body = r.min()
    D = 2 * r_body if r_body > 0.01 else 1.0

    # Body center (approximate as origin for standard cases)
    cx = x[r < r_body * 1.5].mean() if np.any(r < r_body * 1.5) else 0.0
    cy = y[r < r_body * 1.5].mean() if np.any(r < r_body * 1.5) else 0.0

    # View limits clipped to mesh extent
    xlim = (max(cx - 4 * D, x.min()), min(cx + 40 * D, x.max()))
    ylim = (max(cy - 6 * D, y.min()), min(cy + 6 * D, y.max()))

    return xlim, ylim


# ---------------------------------------------------------------------------
# Interpolation grid (built once, shared by all panels)
# ---------------------------------------------------------------------------

class Grid:
    """Uniform interpolation grid for seamless SEM visualization."""

    def __init__(self, x, y, xlim, ylim, ngrid=400):
        aspect = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])
        nxi = ngrid
        nyi = max(int(ngrid * aspect), 50)
        self.xi = np.linspace(xlim[0], xlim[1], nxi)
        self.yi = np.linspace(ylim[0], ylim[1], nyi)
        self.Xi, self.Yi = np.meshgrid(self.xi, self.yi)
        self._pts = np.column_stack([x, y])

        # Pre-compute body mask (auto-detect cylinder radius from mesh)
        r_min = np.sqrt(x**2 + y**2).min()
        self._body_r = r_min * 1.1 if r_min > 0.01 else 0

    def interp(self, vals):
        """Interpolate scattered data to uniform grid, masked inside body."""
        Z = griddata(self._pts, vals, (self.Xi, self.Yi),
                     method='cubic', fill_value=0.0)
        if self._body_r > 0:
            r = np.sqrt(self.Xi**2 + self.Yi**2)
            Z = np.ma.masked_where(r < self._body_r, Z)
        return Z


# ---------------------------------------------------------------------------
# Generic mode figure
# ---------------------------------------------------------------------------

def plot_mode_figure(grid, modes, output_name, suptitle=None):
    """Plot an array of mode panels.

    Parameters
    ----------
    grid   : Grid object
    modes  : list of rows, each row is a list of (Z_interpolated, title) tuples
    output_name : filename (saved in script directory)
    suptitle : optional figure title
    """
    nrows = len(modes)
    ncols = max(len(row) for row in modes)

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(6 * ncols, 2.5 * nrows + 0.5),
                             squeeze=False)

    theta = np.linspace(0, 2 * np.pi, 100)
    cyl_x = grid._body_r / 1.1 * np.cos(theta) if grid._body_r > 0 else None

    for i, row in enumerate(modes):
        # Shared vmax across columns in this row
        vals = [np.abs(Z.compressed()) for Z, _ in row if Z.compressed().size > 0]
        vmax = max((np.percentile(v, 99) for v in vals), default=1.0)
        if vmax < 1e-12:
            vmax = 1.0

        for j, (Z, title) in enumerate(row):
            ax = axes[i, j]
            pcm = ax.pcolormesh(grid.xi, grid.yi, Z, shading='gouraud',
                                cmap=CMAP, vmin=-vmax, vmax=vmax)
            if cyl_x is not None:
                r_body = grid._body_r / 1.1
                ax.fill(r_body * np.cos(theta), r_body * np.sin(theta),
                        'white', ec='k', lw=0.8, zorder=10)
            ax.set_aspect('equal')
            ax.set_xlim(grid.xi[0], grid.xi[-1])
            ax.set_ylim(grid.yi[0], grid.yi[-1])
            ax.set_title(title, fontsize=9)
            fig.colorbar(pcm, ax=ax, shrink=0.8, pad=0.02)

            if j == 0:
                ax.set_ylabel('y')
            if i == nrows - 1:
                ax.set_xlabel('x')

        # Hide unused columns
        for j in range(len(row), ncols):
            axes[i, j].set_visible(False)

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, y=1.01)

    fig.tight_layout()
    out = HERE / output_name
    fig.savefig(out, dpi=DPI, bbox_inches='tight')
    print(f'  Saved {out.name}')
    plt.close(fig)


def plot_mean_figure(grid, x, y, ux, uy):
    """Plot mean flow (velocity magnitude, sequential colormap)."""
    vel = np.sqrt(ux**2 + uy**2)
    Z = grid.interp(vel)

    fig, ax = plt.subplots(figsize=(12, 3.5))
    vmax = min(1.5, np.percentile(np.abs(Z.compressed()), 99.5))
    pcm = ax.pcolormesh(grid.xi, grid.yi, Z, shading='gouraud',
                        cmap=CMAP_SEQ, vmin=0, vmax=vmax)
    if grid._body_r > 0:
        theta = np.linspace(0, 2 * np.pi, 100)
        r = grid._body_r / 1.1
        ax.fill(r * np.cos(theta), r * np.sin(theta),
                'white', ec='k', lw=0.8, zorder=10)
    fig.colorbar(pcm, ax=ax, shrink=0.8, pad=0.02, label=r'$|\mathbf{u}|$')
    ax.set_xlim(grid.xi[0], grid.xi[-1])
    ax.set_ylim(grid.yi[0], grid.yi[-1])
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Mean flow')
    fig.tight_layout()
    out = HERE / 'modes_mean.png'
    fig.savefig(out, dpi=DPI, bbox_inches='tight')
    print(f'  Saved {out.name}')
    plt.close(fig)


# ---------------------------------------------------------------------------
# Decomposition-specific logic
# ---------------------------------------------------------------------------

def do_mean(grid, x, y):
    files = find_mode_files('mea')
    if not files:
        print('  No mean flow files found')
        return
    print(f'  Mean: {files[0].name}')
    _, _, ux, uy = extract_raw(readnek(str(files[0])))
    plot_mean_figure(grid, x, y, ux, uy)


def do_pod(grid, x, y, nmodes):
    files = find_mode_files('pod')
    if not files:
        print('  No POD mode files found')
        return
    n = min(nmodes, len(files))
    print(f'  POD: {len(files)} available, plotting {n}')

    modes = []
    for i in range(n):
        _, _, ux, uy = extract_raw(readnek(str(files[i])))
        Zx = grid.interp(ux)
        Zy = grid.interp(uy)
        modes.append([(Zx, f'POD mode {i+1} — $u_x$'),
                       (Zy, f'POD mode {i+1} — $u_y$')])

    plot_mode_figure(grid, modes, 'modes_pod.png')


def do_dmd(grid, x, y, nmodes):
    files_re = find_mode_files('dm1')
    files_im = find_mode_files('dm2')
    if not files_re:
        print('  No DMD mode files found')
        return
    n = min(nmodes, len(files_re))
    has_im = len(files_im) >= n
    print(f'  DMD: {len(files_re)} available, plotting {n}')

    modes = []
    for i in range(n):
        _, _, _, uy_re = extract_raw(readnek(str(files_re[i])))
        Zy_re = grid.interp(uy_re)
        row = [(Zy_re, f'DMD mode {i+1} — Re($u_y$)')]

        if has_im:
            _, _, _, uy_im = extract_raw(readnek(str(files_im[i])))
            Zy_im = grid.interp(uy_im)
            row.append((Zy_im, f'DMD mode {i+1} — Im($u_y$)'))

        modes.append(row)

    plot_mode_figure(grid, modes, 'modes_dmd.png')


def do_spod(grid, x, y, nmodes):
    files_re = find_mode_files('sRe')
    files_im = find_mode_files('sIm')
    if not files_re:
        print('  No SPOD mode files found')
        return

    nfreqs = len(files_re)
    has_im = len(files_im) >= nfreqs

    # Read spectrum to rank frequencies by energy
    spec_order = list(range(nfreqs))  # fallback: file order
    St_values = None

    for spec_file in ['spod_stream_spectrum.dat', 'spod_spectrum.dat']:
        if (HERE / spec_file).exists():
            data = np.genfromtxt(HERE / spec_file, comments='#')
            if data.size > 0:
                St_values = data[:, 0]
                eig_leading = data[:, 1]
                # Sort by leading eigenvalue, descending (most energetic first)
                spec_order = np.argsort(eig_leading)[::-1].tolist()
            break

    n = min(nmodes, nfreqs)
    print(f'  SPOD: {nfreqs} frequencies, plotting top {n} by energy')

    modes = []
    for rank in range(n):
        fi = spec_order[rank]  # file index (0-based)
        if fi >= nfreqs:
            continue

        freq_label = f'St = {St_values[fi]:.3f}' if St_values is not None else f'freq {fi+1}'

        _, _, _, uy_re = extract_raw(readnek(str(files_re[fi])))
        Zy_re = grid.interp(uy_re)
        row = [(Zy_re, f'SPOD mode {rank+1} ({freq_label}) — Re($u_y$)')]

        if has_im:
            _, _, _, uy_im = extract_raw(readnek(str(files_im[fi])))
            Zy_im = grid.interp(uy_im)
            # Only show Im panel if non-zero (DC/Nyquist are purely real)
            if np.abs(Zy_im.compressed()).max() > 1e-10:
                row.append((Zy_im, f'SPOD mode {rank+1} ({freq_label}) — Im($u_y$)'))

        modes.append(row)

    plot_mode_figure(grid, modes, 'modes_spod.png')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    args = sys.argv[1:]
    flags = set()
    nmodes = 2
    xlim = None
    ylim = None

    i = 0
    while i < len(args):
        a = args[i]
        if   a == '--pod':   flags.add('pod')
        elif a == '--dmd':   flags.add('dmd')
        elif a == '--spod':  flags.add('spod')
        elif a == '--mean':  flags.add('mean')
        elif a == '-n' and i + 1 < len(args):
            i += 1; nmodes = int(args[i])
        elif a == '--xlim' and i + 1 < len(args):
            i += 1; x0, x1 = map(float, args[i].split(',')); xlim = (x0, x1)
        elif a == '--ylim' and i + 1 < len(args):
            i += 1; y0, y1 = map(float, args[i].split(',')); ylim = (y0, y1)
        i += 1

    if not flags:
        flags = {'mean', 'pod', 'dmd', 'spod'}

    return flags, nmodes, xlim, ylim


def main():
    flags, nmodes, xlim, ylim = parse_args()
    print(f'Modal analysis — spatial modes (n={nmodes})')
    print('=' * 50)

    # Read mesh once from any available mode file
    for prefix in ['mea', 'pod', 'dm1', 'sRe']:
        files = find_mode_files(prefix)
        if files:
            field = readnek(str(files[0]))
            x, y, _, _ = extract_raw(field)
            break
    else:
        print('No mode files found.')
        sys.exit(1)

    # Auto-detect limits if not specified
    if xlim is None or ylim is None:
        xl_auto, yl_auto = auto_limits(x, y)
        if xlim is None:
            xlim = xl_auto
        if ylim is None:
            ylim = yl_auto

    # Build shared interpolation grid
    grid = Grid(x, y, xlim, ylim)

    if 'mean' in flags:  do_mean(grid, x, y)
    if 'pod'  in flags:  do_pod(grid, x, y, nmodes)
    if 'dmd'  in flags:  do_dmd(grid, x, y, nmodes)
    if 'spod' in flags:  do_spod(grid, x, y, nmodes)

    print('=' * 50)


if __name__ == '__main__':
    main()
