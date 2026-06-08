"""
spectrum.py — parsers for nekStab eigenvalue / spectrum output files.

Supported formats
-----------------
Generic 2-col  (real, imag)
Spectre_*.dat  4-col: (real, imag, magnitude, flag)  — transient-growth output
dmd_spectrum   8-col with '# DMD' header: mu_real in col5, mu_imag in col6
pod_spectrum   variable cols with '# POD' header: eigenvalue in col1, imag = 0

Public API
----------
parse_spectrum_file(path)           -> np.ndarray, shape (N, 2)
detect_unit_circle_modes(eigvals)   -> np.ndarray of indices
leading_modes(eigvals, k=5)         -> np.ndarray of indices
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _read_header_and_rows(
    path: str,
) -> tuple[str, np.ndarray]:
    """Return (concatenated_header_lines, numeric_data_array).

    Header lines are those starting with '#'.  Subsequent non-blank lines
    must be numeric; non-numeric lines after the first numeric row are
    silently skipped (handles trailing comment blocks).
    """
    header_parts: list[str] = []
    rows: list[list[float]] = []
    past_header = False

    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                header_parts.append(line)
                continue
            past_header = True  # noqa: F841 — kept for clarity
            try:
                tokens = [float(t) for t in line.split()]
                rows.append(tokens)
            except ValueError:
                continue  # non-numeric line after header

    header = " ".join(header_parts).upper()
    if not rows:
        raise ValueError(f"No numeric data found in {path!r}")
    data = np.array(rows, dtype=float)
    return header, data


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_spectrum_file(path: str) -> np.ndarray:
    """Parse an eigenvalue / spectrum file and return complex eigenvalues.

    Returns
    -------
    np.ndarray, shape (N, 2)
        Column 0 = real part, column 1 = imaginary part.

    Format detection
    ----------------
    1. Header contains 'DMD'  → 8-col DMD format; use col5 (mu_real) and
       col6 (mu_imag).
    2. Header contains 'POD'  → variable-col POD format; use col1 as real,
       imag = 0.
    3. ncols == 4              → Spectre_*.dat (real, imag, |mu|, flag);
       use col0 and col1.
    4. ncols == 2              → generic (real, imag); use col0 and col1.
    5. Fallback                → col0 as real, imag = 0.
    """
    header, data = _read_header_and_rows(path)
    ncols = data.shape[1]

    if "DMD" in header:
        # dmd_spectrum.dat: mode |mu| sigma omega St mu_real mu_imag ||Phi||
        real = data[:, 5]
        imag = data[:, 6]

    elif "POD" in header or "SPOD" in header:
        # pod_spectrum.dat: mode eigenvalue energy% cumulative% ||Phi||
        real = data[:, 1]
        imag = np.zeros(len(real))

    elif ncols >= 4:
        # Spectre_Hp.dat / Spectre_NSp.dat: real imag magnitude flag
        real = data[:, 0]
        imag = data[:, 1]

    elif ncols == 2:
        real = data[:, 0]
        imag = data[:, 1]

    else:
        # single-column or unknown
        real = data[:, 0]
        imag = np.zeros(len(real))

    return np.column_stack([real, imag])


def detect_unit_circle_modes(
    eigvals: np.ndarray,
    tol: float = 1e-3,
) -> np.ndarray:
    """Return indices of eigenvalues within *tol* of the unit circle.

    Parameters
    ----------
    eigvals : shape (N, 2) array of [real, imag] pairs.
    tol     : absolute tolerance on ``| |eigval| - 1 |`` (default 1e-3).

    Returns
    -------
    1-D integer array of row indices in *eigvals*.
    """
    eigvals = np.asarray(eigvals, dtype=float)
    magnitudes = np.hypot(eigvals[:, 0], eigvals[:, 1])
    return np.where(np.abs(magnitudes - 1.0) <= tol)[0]


def leading_modes(eigvals: np.ndarray, k: int = 5) -> np.ndarray:
    """Return indices of the top-*k* eigenvalues by magnitude, descending.

    Parameters
    ----------
    eigvals : shape (N, 2) array of [real, imag] pairs.
    k       : number of leading modes to return (clamped to N if k > N).

    Returns
    -------
    1-D integer array of length min(k, N), ordered largest → smallest.
    """
    eigvals = np.asarray(eigvals, dtype=float)
    magnitudes = np.hypot(eigvals[:, 0], eigvals[:, 1])
    k = min(k, len(magnitudes))
    # argsort ascending → reverse for descending
    return np.argsort(magnitudes)[::-1][:k]
