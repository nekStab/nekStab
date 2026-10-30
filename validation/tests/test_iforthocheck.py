"""Orthogonality metric of a Krylov basis (src/ortho.f90).

An orthonormal basis gives a mean absolute inner product below 1e-14.
One repeated vector gives 1/(k-1). The eigensolver check uses the
vector count of the last Arnoldi pass, so a stale later column is
left out. Needs gfortran. The tests skip when it is missing.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"


@pytest.fixture(scope="module")
def driver(build_driver) -> Path:
    return build_driver(
        "ortho",
        [ROOT / "src" / "ortho.f90", FORTRAN / "ortho_driver.f90"],
    )


def ask(exe: Path, text: str) -> list[str]:
    run = subprocess.run([str(exe)], input=text, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    return run.stdout.split()


def basis_text(vectors: np.ndarray, kk: int) -> str:
    """vectors shape (m, k_dim), columns are the Krylov vectors."""
    m, k_dim = vectors.shape
    flat = " ".join(repr(float(x)) for x in vectors.ravel(order="F"))
    return f"basis\n{m} {k_dim} {kk}\n{flat}\n"


def test_orthonormal_basis_is_below_roundoff(driver):
    eye = np.eye(4)
    col, metric = ask(driver, basis_text(eye, 4))
    assert int(col) == 4
    assert float(metric) < 1e-14

    rng = np.random.default_rng(7)
    q, _ = np.linalg.qr(rng.normal(size=(8, 5)))
    col, metric = ask(driver, basis_text(q, 5))
    assert int(col) == 5
    assert float(metric) < 1e-14


def test_one_repeated_vector_gives_one_over_k_minus_one(driver):
    # Q(4) = Q(2). The other products are zero, so the mean is 1/3.
    dots = "0 1 0"
    (metric,) = ask(driver, f"inner\n3\n{dots}\n")
    assert abs(float(metric) - 1.0 / 3.0) < 1e-14

    basis = np.eye(4)
    basis[:, 3] = basis[:, 1]
    col, metric = ask(driver, basis_text(basis, 4))
    assert int(col) == 4
    assert abs(float(metric) - 1.0 / 3.0) < 1e-14


def test_last_pass_ignores_a_stale_column(driver):
    # Columns 1..3 are orthonormal. Column 4 repeats column 1 and belongs
    # to an earlier pass. The check of pass length 3 does not see it.
    basis = np.eye(4)
    basis[:, 3] = basis[:, 0]
    col, metric = ask(driver, basis_text(basis, 3))
    assert int(col) == 3
    assert float(metric) < 1e-14

    for kk in (1, 2, 7, 10):
        (n,) = ask(driver, f"length\n{kk}\n")
        assert int(n) == kk
