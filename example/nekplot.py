#!/usr/bin/env python3
"""nekplot.py: Shared visualization utilities for nekStab 2D example cases.

Provides reusable functions for:
- Reading Nek5000 field files (pymech wrapper)
- Interpolating SEM data to uniform grids
- Detecting and masking solid bodies
- Plotting scalar fields, velocity magnitudes, eigenvalue spectra
- Plotting residual convergence histories
- Plotting OTD Lyapunov exponents
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import griddata
from scipy.spatial import Delaunay


# ── Style ──────────────────────────────────────────────────────────────

def configure_style():
    """Set consistent matplotlib rcParams for all case scripts."""
    plt.rcParams.update({
        'text.usetex': False,
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'legend.fontsize': 7,
        'legend.handlelength': 1.5,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'figure.dpi': 150,
        'savefig.dpi': 600,
        'savefig.format': 'png',
    })


# ── Field I/O ──────────────────────────────────────────────────────────

def read_field(fname):
    """Read a Nek5000 field file via pymech.

    Returns (x, y, fields) where x, y are 1-D coordinate arrays and
    fields is a dict with keys 'vx', 'vy', 'vz' (if present), 'p',
    't' (temperature, if present).
    """
    from pymech.neksuite import readnek

    fname = str(fname)
    field = readnek(fname)
    nel = field.nel
    nx, ny, nz = field.lr1
    npts = nx * ny * max(nz, 1)
    ntotal = nel * npts

    x = np.empty(ntotal)
    y = np.empty(ntotal)
    fields = {}

    for ie, elem in enumerate(field.elem):
        s = slice(ie * npts, (ie + 1) * npts)
        x[s] = elem.pos[0].ravel()
        y[s] = elem.pos[1].ravel()

    # Velocity
    has_vel = any(np.any(elem.vel[0]) for elem in field.elem[:min(5, nel)])
    if has_vel:
        for comp, key in enumerate(['vx', 'vy', 'vz']):
            if comp >= len(field.elem[0].vel):
                break
            arr = np.empty(ntotal)
            for ie, elem in enumerate(field.elem):
                s = slice(ie * npts, (ie + 1) * npts)
                arr[s] = elem.vel[comp].ravel()
            fields[key] = arr

    # Pressure
    has_pres = any(len(elem.pres) > 0 and np.any(elem.pres[0])
                   for elem in field.elem[:min(5, nel)])
    if has_pres:
        arr = np.empty(ntotal)
        for ie, elem in enumerate(field.elem):
            s = slice(ie * npts, (ie + 1) * npts)
            arr[s] = elem.pres[0].ravel()
        fields['p'] = arr

    # Temperature / passive scalars
    if hasattr(field.elem[0], 'temp') and len(field.elem[0].temp) > 0:
        nscalars = len(field.elem[0].temp)
        for iscalar in range(nscalars):
            has_s = any(np.any(elem.temp[iscalar])
                        for elem in field.elem[:min(5, nel)])
            if has_s:
                arr = np.empty(ntotal)
                for ie, elem in enumerate(field.elem):
                    s = slice(ie * npts, (ie + 1) * npts)
                    arr[s] = elem.temp[iscalar].ravel()
                if iscalar == 0:
                    fields['t'] = arr
                else:
                    fields[f's{iscalar:02d}'] = arr

    return x, y, fields, field.time


def find_fields(pattern, directory='.'):
    """Find field files matching a glob pattern, sorted by name."""
    d = Path(directory)
    files = sorted(d.glob(pattern))
    return files


# ── Interpolation ──────────────────────────────────────────────────────

def sem_to_grid(x, y, fields, xlim, ylim, ngrid=400):
    """Interpolate scattered SEM data to a uniform grid.

    Uses cubic interpolation with linear/nearest fallback for boundary
    regions where cubic produces NaN (e.g., near solid bodies).

    Parameters
    ----------
    x, y : 1-D arrays of SEM nodal coordinates
    fields : dict of {name: 1-D values}
    xlim, ylim : (min, max) tuples for grid extent
    ngrid : number of grid points along x-axis

    Returns
    -------
    xi, yi : 1-D grid vectors
    Xi, Yi : 2-D meshgrid arrays
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


# ── Body detection & masking ───────────────────────────────────────────

def detect_body(x, y, tol=0.05):
    """Detect solid bodies as convex hulls around mesh gaps.

    Returns a list of (cx, cy, hull_points) tuples, where cx, cy are
    the centroid and hull_points is an (N, 2) array of boundary vertices.
    Returns empty list if no bodies detected.
    """
    from scipy.spatial import ConvexHull

    # Grid-based gap detection
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    Lx = xmax - xmin
    Ly = ymax - ymin
    ng = 200
    gx = np.linspace(xmin, xmax, ng)
    gy = np.linspace(ymin, ymax, ng)
    dx = Lx / ng
    dy = Ly / ng

    # Count points in each grid cell
    ix = np.clip(((x - xmin) / Lx * (ng - 1)).astype(int), 0, ng - 1)
    iy = np.clip(((y - ymin) / Ly * (ng - 1)).astype(int), 0, ng - 1)
    occ = np.zeros((ng, ng), dtype=bool)
    occ[ix, iy] = True

    # Find empty cells not on domain boundary (margin of 5 cells)
    margin = 5
    empty = ~occ
    empty[:margin, :] = False
    empty[-margin:, :] = False
    empty[:, :margin] = False
    empty[:, -margin:] = False

    if not empty.any():
        return []

    # Connected-component labeling of empty cells
    from scipy.ndimage import label
    labeled, n_bodies = label(empty)

    domain_area = Lx * Ly
    bodies = []
    for b in range(1, n_bodies + 1):
        cells = np.argwhere(labeled == b)
        if len(cells) < 4:
            continue
        # Convert cell indices back to physical coordinates
        bx = gx[cells[:, 0]]
        by = gy[cells[:, 1]]
        cx, cy = bx.mean(), by.mean()
        pts2d = np.column_stack([bx, by])
        try:
            hull = ConvexHull(pts2d)
            # Filter: body area must be < 10% of domain area
            if hull.volume > 0.1 * domain_area:
                continue
            hull_pts = pts2d[hull.vertices]
            bodies.append((cx, cy, hull_pts))
        except Exception:
            pass

    return bodies


def mask_body(Xi, Yi, bodies, margin=0.02):
    """Create a boolean mask for grid points inside detected bodies."""
    if not bodies:
        return np.zeros(Xi.shape, dtype=bool)

    from matplotlib.path import Path as MplPath

    mask = np.zeros(Xi.shape, dtype=bool)
    pts = np.column_stack([Xi.ravel(), Yi.ravel()])

    for _cx, _cy, hull_pts in bodies:
        # Expand hull slightly for cleaner masking
        cx, cy = hull_pts.mean(axis=0)
        expanded = cx + (1 + margin) * (hull_pts[:, 0] - cx), \
                   cy + (1 + margin) * (hull_pts[:, 1] - cy)
        poly = np.column_stack(expanded)
        path = MplPath(poly)
        inside = path.contains_points(pts).reshape(Xi.shape)
        mask |= inside

    return mask


# ── Field plotting ─────────────────────────────────────────────────────

def plot_scalar(ax, xi, yi, Z, bodies=None, cmap='RdBu_r',
                symmetric=True, label=None, vmin=None, vmax=None):
    """Plot a scalar field with pcolormesh and colorbar.

    Parameters
    ----------
    symmetric : if True, use ±vmax for diverging colormap (modes)
                if False, use data range (magnitudes)
    """
    Xi, Yi = np.meshgrid(xi, yi)
    Zm = Z.copy()
    if bodies:
        m = mask_body(Xi, Yi, bodies)
        Zm = np.ma.masked_where(m, Zm)

    if vmin is None or vmax is None:
        p99 = np.nanpercentile(np.abs(Z), 99)
        if symmetric:
            vmin, vmax = -p99, p99
        else:
            vmin, vmax = 0, p99

    pcm = ax.pcolormesh(xi, yi, Zm, shading='gouraud', cmap=cmap,
                        vmin=vmin, vmax=vmax, rasterized=True)
    cb = plt.colorbar(pcm, ax=ax, shrink=0.8, pad=0.02)
    if label:
        cb.set_label(label)

    if bodies:
        add_body_patch(ax, bodies)

    ax.set_xlim(xi[0], xi[-1])
    ax.set_ylim(yi[0], yi[-1])
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')


def plot_velocity_mag(ax, xi, yi, Ux, Uy, bodies=None,
                      cmap='Blues', label=r'$|\mathbf{u}|$'):
    """Plot velocity magnitude sqrt(Ux^2 + Uy^2)."""
    vel_mag = np.sqrt(Ux**2 + Uy**2)
    vmax = min(np.nanpercentile(vel_mag, 99.5),
               np.nanmax(vel_mag))
    plot_scalar(ax, xi, yi, vel_mag, bodies=bodies, cmap=cmap,
                symmetric=False, label=label, vmin=0, vmax=vmax)


def add_streamlines(ax, xi, yi, Ux, Uy, bodies=None, density=1.0):
    """Overlay streamlines on an existing axes."""
    Xi, Yi = np.meshgrid(xi, yi)
    Ux_s, Uy_s = Ux.copy(), Uy.copy()
    if bodies:
        m = mask_body(Xi, Yi, bodies, margin=0.05)
        Ux_s = np.ma.masked_where(m, Ux_s)
        Uy_s = np.ma.masked_where(m, Uy_s)
    try:
        ax.streamplot(xi, yi, Ux_s, Uy_s, color='k', linewidth=0.4,
                      density=density, arrowsize=0.5, arrowstyle='->')
    except Exception:
        pass  # streamplot can fail on masked data


def add_body_patch(ax, bodies):
    """Fill body interiors white with black edge."""
    for _cx, _cy, hull_pts in bodies:
        from matplotlib.patches import Polygon
        poly = Polygon(hull_pts, closed=True, fc='white', ec='k',
                       lw=0.8, zorder=10)
        ax.add_patch(poly)


# ── Spectrum plotting ──────────────────────────────────────────────────

def _read_spectrum(dat_file):
    """Read a Spectre_*.dat file. Returns (col0, col1, residual, flag).

    Format: 3-4 columns of Fortran-formatted floats.
    Col 0: sigma (NS) or Re(mu) (H)
    Col 1: omega (NS) or Im(mu) (H)
    Col 2: residual (optional)
    Col 3: converged flag (optional, 1=converged)
    """
    data = np.genfromtxt(str(dat_file))
    if data.ndim == 1:
        data = data.reshape(1, -1)
    c0 = data[:, 0]
    c1 = data[:, 1]
    res = data[:, 2] if data.shape[1] > 2 else np.zeros(len(c0))
    flag = data[:, 3] if data.shape[1] > 3 else np.ones(len(c0))
    return c0, c1, res, flag


def plot_spectrum_H(ax, dat_file, tolerance=1e-6):
    """Plot eigenvalues on the unit circle (Floquet/Hessenberg).

    Converged eigenvalues (residual < tolerance) in black,
    unconverged in gray.
    """
    dat_file = Path(dat_file)
    if not dat_file.exists():
        ax.text(0.5, 0.5, f'{dat_file.name}\nnot found',
                transform=ax.transAxes, ha='center', va='center')
        return

    Re_mu, Im_mu, res, _flag = _read_spectrum(dat_file)

    # Unit circle
    theta = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(theta), np.sin(theta), 'r-', lw=0.5, zorder=1)

    conv = res < tolerance
    if conv.any():
        ax.scatter(Re_mu[conv], Im_mu[conv], s=12, marker='o',
                   facecolors='k', edgecolors='k', linewidths=0.3,
                   zorder=3, label='converged')
    if (~conv).any():
        ax.scatter(Re_mu[~conv], Im_mu[~conv], s=8, marker='o',
                   facecolors='lightgray', edgecolors='gray',
                   linewidths=0.3, alpha=0.5, zorder=2, label='unconverged')

    ax.axhline(0, lw=0.3, color='k', ls=':')
    ax.axvline(0, lw=0.3, color='k', ls=':')
    ax.set_xlabel(r'$\Re(\mu)$')
    ax.set_ylabel(r'$\Im(\mu)$')
    ax.set_aspect('equal')
    ax.legend(fontsize=6, loc='best')


def plot_spectrum_NS(ax, dat_file, tolerance=1e-6, freq=True):
    """Plot eigenvalues as sigma vs omega (or frequency f=omega/2pi).

    Converged eigenvalues in black, unconverged in gray.
    """
    dat_file = Path(dat_file)
    if not dat_file.exists():
        ax.text(0.5, 0.5, f'{dat_file.name}\nnot found',
                transform=ax.transAxes, ha='center', va='center')
        return

    sigma, omega, res, _flag = _read_spectrum(dat_file)
    x = omega / (2 * np.pi) if freq else omega
    xlabel = r'$f = \omega/2\pi$' if freq else r'$\omega$'

    conv = res < tolerance
    if conv.any():
        ax.scatter(x[conv], sigma[conv], s=12, marker='o',
                   facecolors='k', edgecolors='k', linewidths=0.3,
                   zorder=3, label='converged')
    if (~conv).any():
        ax.scatter(x[~conv], sigma[~conv], s=8, marker='o',
                   facecolors='lightgray', edgecolors='gray',
                   linewidths=0.3, alpha=0.5, zorder=2, label='unconverged')

    ax.axhline(0, lw=0.3, color='k', ls=':')
    ax.axvline(0, lw=0.3, color='k', ls=':')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$\sigma$')
    ax.legend(fontsize=6, loc='best')


# ── Convergence plotting ──────────────────────────────────────────────

def plot_residuals(ax, dat_file, label=None, color='k'):
    """Plot SFD/BoostConv residuals (residu.dat: time, L2, L2_t)."""
    dat_file = Path(dat_file)
    if not dat_file.exists():
        ax.text(0.5, 0.5, f'{dat_file.name}\nnot found',
                transform=ax.transAxes, ha='center', va='center')
        return

    data = np.genfromtxt(str(dat_file))
    if data.ndim == 1:
        data = data.reshape(1, -1)
    t = data[:, 0]
    r = data[:, 1]

    ax.semilogy(t, r, color=color, lw=0.8, label=label)
    ax.axhline(r.min(), color=color, lw=0.4, ls='--', alpha=0.5)
    ax.set_xlabel('t')
    ax.set_ylabel(r'$\|r\|$')
    if label:
        ax.legend(fontsize=6)


def plot_newton_convergence(ax, directory='.'):
    """Plot Newton + Arnoldi convergence from residu_*.dat files."""
    d = Path(directory)

    # Arnoldi residuals (continuous curve)
    arnoldi_file = d / 'residu_arnoldi.dat'
    if arnoldi_file.exists():
        data = np.genfromtxt(str(arnoldi_file))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        # Columns: k(0), time(1), tol(2), beta2(3), dtol(4)
        k = np.arange(len(data))
        res_arn = data[:, 3]   # beta2 = Arnoldi residual
        tol_arn = data[:, 2]   # solver tolerance
        ax.semilogy(k, res_arn, 'k-', lw=0.6, label='Arnoldi', zorder=3)
        ax.semilogy(k, tol_arn, 'b--', lw=0.5, label='tol', zorder=2)

    # Newton residuals (markers at k_sum positions)
    newton_file = d / 'residu_newton.dat'
    if newton_file.exists():
        data = np.genfromtxt(str(newton_file))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        # Columns: i, total_calls, iter_calls, k_sum, tottime, solver_tol, residual, dtol
        k_sum = data[:, 3]
        res_nwt = data[:, 6]
        dtol = data[0, 7]
        ax.semilogy(k_sum, res_nwt, 'bs', ms=4, mfc='none', lw=0.8,
                    label='Newton', zorder=4)
        ax.axhline(dtol, color='r', ls='--', lw=0.5, label='dtol', zorder=1)

    ax.set_xlabel('Arnoldi iterations')
    ax.set_ylabel(r'$\|r\|^2$')
    ax.legend(fontsize=6, ncol=2)
    ax.grid(True, which='both', ls=':', lw=0.3, alpha=0.5)


def plot_tdf_residuals(ax, dat_file, label=None, color='k'):
    """Plot TDF residuals (residu_tdf.dat: same format as residu.dat)."""
    plot_residuals(ax, dat_file, label=label, color=color)


# ── OTD plotting ───────────────────────────────────────────────────────

def plot_otd_exponents(ax, directory='.'):
    """Plot OTD Lyapunov exponents from otd_ftle/eigenvalues/growth_rates.dat files."""
    d = Path(directory)
    colors_cycle = ['b', 'g', 'r', 'c', 'm']

    for fname, ls, prefix, suffix in [
        ('otd_ftle.dat', '-',  r'$\lambda', '$'),
        ('otd_eigenvalues.dat', '--', r'$\Re(\lambda', ')$'),
        ('otd_growth_rates.dat', '-.', r'$\sigma', '$'),
    ]:
        fpath = d / fname
        if not fpath.exists():
            continue
        data = np.genfromtxt(str(fpath))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        t = data[:, 0]
        n_cols = data.shape[1] - 1
        for j in range(n_cols):
            label = prefix + f'_{j+1}' + suffix
            c = colors_cycle[j % len(colors_cycle)]
            ax.plot(t, data[:, j + 1], ls=ls, lw=0.8, color=c, label=label)

    ax.set_xlabel('t')
    ax.set_ylabel(r'$\lambda$, $\Re(\lambda)$, $\sigma$')
    ax.legend(fontsize=5, ncol=3, loc='best')
    ax.grid(True, ls=':', lw=0.3, alpha=0.5)


# ── Transient growth ───────────────────────────────────────────────────

def plot_transient_growth(ax, directory='.', ref_file=None):
    """Plot G(t) = max energy amplification envelope.

    Reads Spectre_Hp.dat from subdirectories t_* to get G(t) = max eigenvalue.
    Optionally overlays a reference curve from ref_file.
    """
    d = Path(directory)

    if ref_file and Path(ref_file).exists():
        ref = np.genfromtxt(str(ref_file))
        if ref.ndim == 1:
            ref = ref.reshape(1, -1)
        ax.plot(ref[:, 0], ref[:, 1], 'k-', lw=0.5, label='Reference')

    # Scan t_* directories for Spectre_Hp.dat
    t_vals = []
    g_vals = []
    for td in sorted(d.glob('t_*')):
        sp = td / 'Spectre_Hp.dat'
        if sp.exists():
            try:
                t_val = float(td.name.split('_')[1])
                data = np.genfromtxt(str(sp))
                if data.ndim == 1:
                    data = data.reshape(1, -1)
                g_vals.append(data[0, 0])  # Leading eigenvalue
                t_vals.append(t_val)
            except Exception:
                continue

    if t_vals:
        ax.scatter(t_vals, g_vals, s=15, marker='o', facecolors='none',
                   edgecolors='k', linewidths=0.5, zorder=3, label='nekStab')

    ax.set_xlabel('t')
    ax.set_ylabel('G(t) = E(t)/E(0)')
    ax.legend(fontsize=6)


# ── Convenience ────────────────────────────────────────────────────────

def load_and_interp(field_file, field_keys, xlim, ylim, ngrid=400):
    """Read a field file and interpolate specified components to a grid.

    Parameters
    ----------
    field_file : path to .f????? file
    field_keys : list of field component names to interpolate (e.g., ['vy'])
    xlim, ylim : plot bounds

    Returns
    -------
    x, y : raw SEM coordinates
    xi, yi, Xi, Yi : grid arrays
    interp : dict of interpolated fields
    bodies : detected body list
    """
    x, y, fields, time = read_field(field_file)
    to_interp = {k: fields[k] for k in field_keys if k in fields}
    if not to_interp:
        return None
    xi, yi, Xi, Yi, interp = sem_to_grid(x, y, to_interp, xlim, ylim, ngrid)
    bodies = detect_body(x, y)
    return x, y, xi, yi, Xi, Yi, interp, bodies, time


# ══════════════════════════════════════════════════════════════════════
# Paper-style plotting (tricontourf, Delaunay, inset colorbars)
# Matches AMR_Krylov_V5/Main/bkp_scripts visual style exactly
# ══════════════════════════════════════════════════════════════════════

# Journal column width: 240.7103 pt = 3.331 inches
COL_WIDTH = 240.7103 / 72.27

def discrete_cmap(N, base_cmap='RdBu'):
    """Create an N-bin discrete colormap from a base colormap (paper style)."""
    base = plt.cm.get_cmap(base_cmap)
    colors = base(np.linspace(0, 1, N))
    return ListedColormap(colors, name=f'{base.name}{N}')


def make_triangulation(x, y):
    """Build Delaunay triangulation from SEM coordinates."""
    tri = Delaunay(np.vstack((x, y)).T)
    return Triangulation(x, y, triangles=tri.simplices)


def tricontourf(ax, triang, q, levels=257, cmap='RdBu', extend='both',
                vmin=None, vmax=None, symmetric=False):
    """Paper-style tricontourf with discrete_cmap.

    Parameters
    ----------
    ax : matplotlib axes
    triang : Triangulation object
    q : 1-D field values at triangle vertices
    levels : int or array of contour levels
    cmap : base colormap name
    extend : 'both', 'max', 'min', 'neither'
    vmin, vmax : explicit bounds (overrides symmetric)
    symmetric : if True, use ±max(|q|) for vmin/vmax

    Returns
    -------
    ContourSet (for colorbar creation)
    """
    if isinstance(levels, int):
        n = levels
        if vmin is not None and vmax is not None:
            lev = np.linspace(vmin, vmax, n)
        elif symmetric:
            bd = np.nanpercentile(np.abs(q), 99)
            lev = np.linspace(-bd, bd, n)
        else:
            lev = np.linspace(np.nanmin(q), np.nanmax(q), n)
    else:
        lev = levels
        n = len(lev)

    cm = discrete_cmap(n, cmap)
    return ax.tricontourf(triang, q, levels=lev, cmap=cm, extend=extend)


def inset_colorbar(ax, mappable, orientation='vertical', width="3%",
                   height="70%", loc=1, ticks=None, tick_labels=None):
    """Add a compact inset colorbar inside the axes (paper style).

    Parameters
    ----------
    ax : parent axes
    mappable : ContourSet or ScalarMappable
    orientation : 'vertical' or 'horizontal'
    width, height : size of inset (swapped for horizontal)
    loc : location code (1=upper right, 2=upper left, etc.)
    ticks : tick values
    tick_labels : custom tick labels (strings)
    """
    # Map convenience names to valid matplotlib orientation
    cb_orient = 'horizontal' if orientation in ('horizontal', 'top', 'bottom') \
        else 'vertical'
    cbaxes = inset_axes(ax, width=width, height=height, loc=loc)
    cbar = plt.colorbar(mappable, cax=cbaxes, orientation=cb_orient,
                        ticks=ticks)
    cbar.ax.tick_params(labelsize=5, length=2, pad=1)
    if cb_orient == 'vertical':
        cbar.ax.get_yaxis().set_ticks_position('left')
    if tick_labels is not None:
        if cb_orient == 'vertical':
            cbar.ax.set_yticklabels(tick_labels, fontsize=5)
        else:
            cbar.ax.set_xticklabels(tick_labels, fontsize=5)
    return cbar


def panel_label(ax, text, x=-0.12, y=1.05):
    """Add a bold panel label like (a), (b) outside the axes."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=9,
            fontweight='bold', va='bottom', ha='right')


# ── Geometry patches ──────────────────────────────────────────────────

def add_annulus_patches(ax, r_inner=1.0, r_outer=2.0):
    """Thermosyphon: inner filled circle + outer ring."""
    inner = mpl.patches.Circle((0, 0), r_inner, facecolor='white',
                                edgecolor='black', lw=0.6, zorder=10)
    outer = mpl.patches.Circle((0, 0), r_outer, fill=False,
                                edgecolor='black', lw=0.6, zorder=10)
    ax.add_artist(outer)
    ax.add_artist(inner)


def add_step_patch(ax, x0=-5, y0=-1, w=5, h=1):
    """Backward-facing step: solid rectangle body."""
    rect = mpl.patches.Rectangle((x0, y0), w, h, facecolor='lightgray',
                                  edgecolor='black', lw=0.3, zorder=10)
    ax.add_artist(rect)


def add_cylinder_patches(ax, cx=0, cy=0, radius=0.5):
    """Single cylinder body patch."""
    patch = mpl.patches.Circle((cx, cy), radius, facecolor='lightgray',
                                edgecolor='black', lw=0.3, zorder=10)
    ax.add_artist(patch)


def add_dual_cylinder_patches(ax, gap=0.7, radius=0.5):
    """Flip-flop: two side-by-side cylinders with gap g."""
    cy = gap / 2.0 + radius
    for sign in (+1, -1):
        patch = mpl.patches.Circle((0, sign * cy), radius,
                                    facecolor='lightgray',
                                    edgecolor='black', lw=0.3, zorder=10)
        ax.add_artist(patch)


def add_naca0012_patch(ax, chord=1.0, cx=0.25, cy=0):
    """NACA 0012 airfoil body patch (symmetric, analytical profile).

    Default cx=0.25 places leading edge at x=-0.25, trailing at x=0.75,
    matching the standard nekStab mesh (quarter-chord at origin).
    """
    t_naca = 0.12  # max thickness ratio
    x_pts = np.linspace(0, 1, 200)
    yt = 5 * t_naca * (0.2969 * np.sqrt(x_pts) - 0.1260 * x_pts
                        - 0.3516 * x_pts**2 + 0.2843 * x_pts**3
                        - 0.1015 * x_pts**4)
    # Close upper + lower surface; shift so that cx is at chord/2
    xs = np.concatenate([x_pts, x_pts[::-1]]) * chord + (cx - chord / 2)
    ys = np.concatenate([yt, -yt[::-1]]) * chord + cy
    poly = mpl.patches.Polygon(np.column_stack([xs, ys]), closed=True,
                                facecolor='lightgray', edgecolor='black',
                                lw=0.3, zorder=10)
    ax.add_artist(poly)


# ── Paper spectrum style ──────────────────────────────────────────────

def plot_spectrum_H_paper(ax, dat_file, tolerance=1e-5, sized=8,
                          symb='o', cor='royalblue', label=None):
    """Paper-style Floquet spectrum: unit disk filled gray, scatter.

    Only plots converged eigenvalues (residual < tolerance).
    """
    dat_file = Path(dat_file)
    if not dat_file.exists():
        ax.text(0.5, 0.5, f'{dat_file.name}\nnot found',
                transform=ax.transAxes, ha='center', va='center')
        return

    Re_mu, Im_mu, res, _flag = _read_spectrum(dat_file)

    # Plot converged eigenvalues only
    first = True
    for k in range(len(Re_mu)):
        if res[k] < tolerance:
            lbl = label if first else None
            ax.scatter(Re_mu[k], Im_mu[k], s=sized, c=cor, marker=symb,
                       edgecolors='k', linewidths=0.3, label=lbl, zorder=3)
            first = False

    ax.axhline(0, lw=0.5, c='gray', ls='dotted')
    ax.axvline(0, lw=0.5, c='gray', ls='dotted')


def setup_unit_circle_axes(ax, lim=1.5):
    """Configure axes for unit-circle Floquet spectrum (paper style)."""
    circle = mpl.patches.Circle((0, 0), 1, color='gray', alpha=0.5,
                                 zorder=0, linewidth=0, fill=True)
    ax.add_artist(circle)
    ax.set_aspect('equal')
    ticks = [-lim, 0, lim]
    tick_labels = [str(t) for t in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(r'$\Re(\mu)$')
    ax.set_ylabel(r'$\Im(\mu)$')


def plot_spectrum_NS_paper(ax, dat_file, tolerance=1e-5, freq=True,
                           sized=14, symb='d', cor='k', label=None):
    """Paper-style NS spectrum: sigma vs f with converged filter."""
    dat_file = Path(dat_file)
    if not dat_file.exists():
        ax.text(0.5, 0.5, f'{dat_file.name}\nnot found',
                transform=ax.transAxes, ha='center', va='center')
        return

    sigma, omega, res, _flag = _read_spectrum(dat_file)
    b = (2.0 * np.pi) if freq else 1.0

    first = True
    for k in range(len(sigma)):
        if res[k] < tolerance:
            lbl = label if first else None
            ax.scatter(omega[k] / b, sigma[k], s=sized, c=cor, marker=symb,
                       edgecolors='k', linewidths=0.3, label=lbl, zorder=3)
            first = False

    ax.axhline(0, lw=0.5, c='gray', ls='dotted')
    ax.axvline(0, lw=0.5, c='gray', ls='dotted')
    ax.set_xlabel(r'$f$')
    ax.set_ylabel(r'$\sigma$')

# ── Reynolds number parsing ────────────────────────────────────────

def re_from_par(par_path):
    """Parse Reynolds number from a Nek5000 .par file or directory.
    
    If par_path is a directory, finds the single *.par file in it.
    Parses the [VELOCITY] section for the 'viscosity' value.
    Handles inline comments (# ...), whitespace, and case-insensitive keys.
    
    Returns Reynolds number as float:
      - If viscosity < 0: Re = -viscosity (negative convention)
      - If viscosity > 0: Re = 1.0 / viscosity (positive convention)
    Returns None if not found or unreadable.
    """
    from pathlib import Path
    
    par_path = Path(par_path)
    
    # If directory, find the single *.par file
    if par_path.is_dir():
        par_files = list(par_path.glob('*.par'))
        if len(par_files) != 1:
            return None
        par_path = par_files[0]
    
    # Try to read and parse
    try:
        with open(par_path, 'r') as f:
            lines = f.readlines()
    except Exception:
        return None
    
    # Find [VELOCITY] section and extract viscosity
    in_velocity = False
    for line in lines:
        # Strip inline comments
        if '#' in line:
            line = line[:line.index('#')]
        line = line.strip()
        
        if line.lower() == '[velocity]':
            in_velocity = True
            continue
        
        # Stop at next section
        if line.startswith('[') and in_velocity:
            break
        
        # Look for viscosity key
        if in_velocity and '=' in line:
            key, val = line.split('=', 1)
            if key.strip().lower() == 'viscosity':
                try:
                    nu = float(val.strip())
                    if nu < 0:
                        return -nu
                    else:
                        return 1.0 / nu if nu != 0 else None
                except ValueError:
                    return None
    
    return None


def re_label(ax, re, loc='upper left'):
    """Annotate axes with Reynolds number label.
    
    If re is not None, adds text "$Re = {re:g}$" to the axes in mathtext.
    Uses fontsize ~7 and places in axes coordinates near top-left by default.
    No-op if re is None.
    
    Parameters
    ----------
    ax : matplotlib axes
    re : float or None
        Reynolds number to display
    loc : str
        Location preset ('upper left', 'upper right', etc.)
    """
    if re is None:
        return
    
    # Map location to axes coordinates
    loc_map = {
        'upper left': (0.02, 0.95),
        'upper right': (0.98, 0.95),
        'lower left': (0.02, 0.05),
        'lower right': (0.98, 0.05),
    }
    xy = loc_map.get(loc, (0.02, 0.95))
    ha = 'left' if xy[0] < 0.5 else 'right'
    va = 'top' if xy[1] > 0.5 else 'bottom'
    
    ax.text(xy[0], xy[1], f'$Re = {re:g}$', transform=ax.transAxes,
            fontsize=7, ha=ha, va=va, bbox=None)
