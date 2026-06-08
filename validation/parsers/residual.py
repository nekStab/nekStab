"""
residual.py — parsers for nekStab / Nek5000 residual and force-history files.

Supported formats
-----------------
2-col  (time, residual)                                    -> boostconv / SFD simple logs
3-col  (time, residual, rate)                              -> residu_sfd.dat / residu_dmt.dat / boostconv
5-col  (step, ?, ?, residual, tol)                         -> residu_arnoldi.dat
7-col  (step, ?, iter, cumul, residual_prev, residual, tol)-> residu_gmres.dat
8-col  (step, cumul, iter, inner, time, res_prev, res, tol)-> residu_newton.dat
11-col (step, time, col1..col8, ...)                       -> lift_drag.dat (Nek5000 force history)

Public API
----------
parse_residual_log(path)   -> (time, residual, rate | None)
parse_force_history(path)  -> dict[str, np.ndarray]
detect_convergence(...)    -> dict
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_numeric_rows(path: str) -> np.ndarray:
    """Read whitespace-delimited numeric rows; skip comment / header lines.

    Lines that start with '#' or contain non-numeric tokens (after substituting
    'Infinity' / 'infinity' with a sentinel) are dropped.
    Returns a 2-D float array, shape (nrows, ncols).  Raises ValueError when
    no valid rows are found.
    """
    rows = []
    with open(path, "r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            # Replace Infinity tokens so float() succeeds
            line = line.replace("Infinity", "inf").replace("infinity", "inf")
            try:
                tokens = [float(t) for t in line.split()]
                rows.append(tokens)
            except ValueError:
                continue  # skip header / label lines
    if not rows:
        raise ValueError(f"No numeric data found in {path!r}")
    return np.array(rows, dtype=float)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_residual_log(
    path: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Parse a nekStab / Nek5000 residual log file.

    Returns
    -------
    time     : 1-D float array (simulation time or iteration index)
    residual : 1-D float array (convergence residual)
    rate     : 1-D float array or None (convergence rate, when present)

    Format dispatch (by column count)
    ----------------------------------
    2 cols  -> (time, residual); rate = None
    3 cols  -> (time, residual, rate); 'Infinity' → np.inf already handled
    5 cols  -> arnoldi: time = col0 (step as float), residual = col3; rate = None
    7 cols  -> gmres: time = col0 (step), residual = col5; rate = None
    8 cols  -> newton: time = col4 (sim time), residual = col6; rate = None
    """
    data = _load_numeric_rows(path)
    ncols = data.shape[1]

    if ncols == 2:
        return data[:, 0], data[:, 1], None

    if ncols == 3:
        rate = data[:, 2].copy()
        return data[:, 0], data[:, 1], rate

    if ncols == 5:
        # residu_arnoldi.dat: step, something, something, residual, tol
        return data[:, 0], data[:, 3], None

    if ncols == 7:
        # residu_gmres.dat: step, ?, iter, cumul, res_prev, residual, tol
        return data[:, 0], data[:, 5], None

    if ncols >= 8:
        # residu_newton.dat: step, cumul, iter, inner, time, res_prev, residual, tol
        return data[:, 4], data[:, 6], None

    # fallback: treat col0 as time, col1 as residual
    return data[:, 0], data[:, 1], None


def parse_force_history(path: str) -> dict[str, np.ndarray]:
    """Parse a Nek5000 force / lift-drag history file.

    Expected format: whitespace-delimited rows where column 0 is the
    integer step counter and column 1 is simulation time; subsequent
    columns are force components.

    Returns
    -------
    dict with keys ``'step'``, ``'time'``, and ``'col_<N>'`` (1-indexed)
    for each force column (N ≥ 2).

    Example (lift_drag.dat, 11 columns total)::

        0  0.0  0.0  0.0  ...  (10 numeric columns after 'step')

    Raises
    ------
    ValueError  if the file contains fewer than 2 numeric columns.
    """
    data = _load_numeric_rows(path)
    if data.shape[1] < 2:
        raise ValueError(
            f"{path!r}: expected ≥ 2 columns (step, time, ...), "
            f"got {data.shape[1]}"
        )
    result: dict[str, np.ndarray] = {
        "step": data[:, 0].astype(int).astype(float),
        "time": data[:, 1],
    }
    for k in range(2, data.shape[1]):
        result[f"col_{k}"] = data[:, k]
    return result


def detect_convergence(
    residual: np.ndarray,
    time: np.ndarray,
    tol: float = 1e-6,
) -> dict:
    """Detect convergence in a residual history.

    Parameters
    ----------
    residual : 1-D array of residual values.
    time     : 1-D array of corresponding times / iteration indices
               (must have the same length as *residual*).
    tol      : convergence threshold (default 1e-6).

    Returns
    -------
    dict with keys:

    ``'converged'``       bool — True if ``residual[-1] <= tol``
    ``'converge_time'``   float or None — first ``time[i]`` where
                          ``residual[i] <= tol``; None if never converged
    ``'final_residual'``  float — ``residual[-1]``
    """
    residual = np.asarray(residual, dtype=float)
    time = np.asarray(time, dtype=float)

    final = float(residual[-1])
    converged = bool(final <= tol)

    converge_time: float | None = None
    idx = np.where(residual <= tol)[0]
    if idx.size > 0:
        converge_time = float(time[idx[0]])

    return {
        "converged": converged,
        "converge_time": converge_time,
        "final_residual": final,
    }
