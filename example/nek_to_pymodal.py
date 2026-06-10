#!/usr/bin/env python3
"""nek_to_pymodal.py: Convert nekStab spectral-element snapshots to pyModal npz format.

Reads a time series of Nek5000 .f field files (unstructured spectral-element mesh),
interpolates velocity fields onto a uniform Cartesian grid, and writes the consolidated
.npz archive in the format expected by pyModal's DNamiXNPZLoader.

The output npz contains:
  - x : 1D array of Cartesian x grid coordinates, shape (Nx,)
  - y : 1D array of Cartesian y grid coordinates, shape (Ny,)
  - u : 3D array of x-velocity snapshots, shape (Ns, Nx, Ny)
  - v : 3D array of y-velocity snapshots, shape (Ns, Nx, Ny)
  - times : 1D array of snapshot times, shape (Ns,)
"""

from pathlib import Path
import sys
import time

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

# Import nekplot utilities
sys.path.insert(0, str(Path(__file__).parent))
import nekplot as nk


def nek_snapshots_to_npz(snapshot_paths, out_path, bbox = (-2.0, 20.0, -4.0, 4.0),
                          nx = 256, ny = 128, fields = ('vx', 'vy'), fill = 0.0):
    """Convert Nek spectral-element snapshots to pyModal npz format.

    Parameters
    ----------
    snapshot_paths : list of Path or str
        Sorted list of Nek .f field file paths (time series).
    out_path : Path or str
        Output .npz file path.
    bbox : tuple of (xmin, xmax, ymin, ymax)
        Bounding box for the Cartesian grid.
    nx, ny : int
        Grid resolution (number of points along x and y).
    fields : tuple of str
        Field component names to extract (e.g., ('vx', 'vy')).
    fill : float
        Fill value for NaN points (e.g., inside solid bodies).

    Returns
    -------
    dict
        Summary dict with keys 'out_path', 'Ns', 'Nx', 'Ny'.
    """

    snapshot_paths = sorted([Path(p) for p in snapshot_paths])
    out_path = Path(out_path)
    Ns = len(snapshot_paths)

    # Build Cartesian grid
    x_grid = np.linspace(bbox[0], bbox[1], nx)
    y_grid = np.linspace(bbox[2], bbox[3], ny)

    # Storage for all snapshots
    u_all = np.zeros((Ns, nx, ny), dtype=np.float64)
    v_all = np.zeros((Ns, nx, ny), dtype=np.float64)
    times = np.zeros(Ns, dtype=np.float64)

    # Reusable triangulation and interpolators (assume fixed mesh)
    tri = None
    interp_u = None
    interp_v = None
    x_prev = None
    y_prev = None

    # Process each snapshot
    for i, snap_path in enumerate(snapshot_paths):
        if (i + 1) % 10 == 0:
            print(f"  Processing snapshot {i + 1}/{Ns}...")

        x, y, snap_fields, snap_time = nk.read_field(snap_path)
        times[i] = snap_time

        # Check if mesh has changed; if so, rebuild triangulation
        if tri is None or len(x) != len(x_prev):
            # Build Delaunay triangulation from scattered nodes
            pts = np.column_stack([x, y])
            tri = Delaunay(pts)
            x_prev = x.copy()
            y_prev = y.copy()
            # Rebuild interpolators for new mesh
            u_vals = snap_fields.get('vx', np.zeros_like(x))
            v_vals = snap_fields.get('vy', np.zeros_like(x))
            interp_u = LinearNDInterpolator(tri, u_vals, fill_value = fill)
            interp_v = LinearNDInterpolator(tri, v_vals, fill_value = fill)
        else:
            # Mesh unchanged: reuse tri, update field values
            u_vals = snap_fields.get('vx', np.zeros_like(x))
            v_vals = snap_fields.get('vy', np.zeros_like(x))
            interp_u = LinearNDInterpolator(tri, u_vals, fill_value = fill)
            interp_v = LinearNDInterpolator(tri, v_vals, fill_value = fill)

        # Create meshgrid for interpolation (indexing='ij' for (nx, ny) shape)
        Xi, Yi = np.meshgrid(x_grid, y_grid, indexing = 'ij')
        pts_grid = np.column_stack([Xi.ravel(), Yi.ravel()])

        # Interpolate onto grid
        u_interp = interp_u(pts_grid).reshape(Xi.shape)
        v_interp = interp_v(pts_grid).reshape(Xi.shape)

        # Replace NaN with fill value
        u_interp = np.nan_to_num(u_interp, nan = fill)
        v_interp = np.nan_to_num(v_interp, nan = fill)

        u_all[i, :, :] = u_interp
        v_all[i, :, :] = v_interp

    # Write npz
    out_path.parent.mkdir(parents = True, exist_ok = True)
    np.savez(out_path, x = x_grid, y = y_grid, u = u_all, v = v_all, times = times)

    return {'out_path': str(out_path), 'Ns': Ns, 'Nx': nx, 'Ny': ny}


def main():
    """Default invocation: convert cylinder_re100 modal POD snapshots."""
    case_dir = Path(__file__).parent / 'cylinder_re100' / '600_modal_pod'
    snapshots = sorted(case_dir.glob('1cyl0.f*'))

    if not snapshots:
        print(f"No snapshots found in {case_dir}")
        return

    out_file = case_dir / 'modal_snapshots.npz'

    print(f"Converting {len(snapshots)} snapshots...")
    print(f"  From: {case_dir}")
    print(f"  To: {out_file}")

    start_time = time.time()
    result = nek_snapshots_to_npz(snapshots, out_file)
    elapsed = time.time() - start_time

    print(f"\nDone in {elapsed:.2f} s")
    print(f"  Ns={result['Ns']}, Nx={result['Nx']}, Ny={result['Ny']}")
    print(f"  Output: {result['out_path']}")


if __name__ == '__main__':
    main()
