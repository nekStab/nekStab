"""Exact zeros of the Krylov weight (src/bm1s_mask.f90).

The mask sets a field to zero where the weight is exactly zero and leaves
every other point unchanged. In the solver the weight is bm1s. The sponge
stores those zeros; the test reads the values it is given.
Needs gfortran; the tests skip when it is missing.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"

# Each case is the weight text, then the field text. A zero, a signed zero,
# and a tiny nonzero sit next to ordinary values so a threshold would fail.
CASES = [
    (
        "0 -0.0 1e-300 1 2.5 0",
        "9 -4 1.25 0 3.5 -8",
    ),
    (
        "1 2 3 4",
        "1.5 -2.5 3.5 -4.5",
    ),
    (
        "0 0 0",
        "1e10 -1e10 7",
    ),
    (
        "1",
        "42",
    ),
    (
        "0",
        "-6.5",
    ),
    (
        "4.5 0 -8 1e-320 0",
        "1 2 3 4 5",
    ),
]


@pytest.fixture(scope="module")
def driver(build_driver) -> Path:
    return build_driver(
        "bm1s_mask",
        [ROOT / "src" / "bm1s_mask.f90", FORTRAN / "bm1s_mask_driver.f90"],
    )


def run_case(exe: Path, weights: str, fields: str) -> tuple[int, list[tuple[int, str, str]]]:
    n = len(weights.split())
    assert n == len(fields.split())
    text = f"{n}\n{weights}\n{fields}\n"
    run = subprocess.run([str(exe)], input=text, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    rows = run.stdout.strip().splitlines()
    assert len(rows) == n + 1
    has_zero = int(rows[0])
    parsed = []
    for row in rows[1:]:
        flag, before, after = row.split()
        parsed.append((int(flag), before.upper(), after.upper()))
    return has_zero, parsed


@pytest.mark.parametrize("weights,fields", CASES)
def test_exact_zeros_are_cleared_and_other_points_stay(driver, weights, fields):
    has_zero, rows = run_case(driver, weights, fields)
    zero = "0000000000000000"
    assert has_zero == int(any(flag == 1 for flag, _, _ in rows))
    for flag, before, after in rows:
        if flag == 1:
            assert after == zero
        else:
            assert after == before
