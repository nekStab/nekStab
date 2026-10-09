"""Ritz value selection of the Krylov-Schur restart (src/ks_select.f90).

The module is compiled alone (it needs no Nek5000) and run once per case by
fortran/ks_select_driver.f90. The cases are built here, not read from a run.

An earlier version removed single values by residual and then completed the
conjugate pairs. The completion returned the count to n, no Arnoldi step
followed, and the residual test of krylov_schur reported n - 2 converged
eigenvalues (an ellipse with k_dim = 160: "keeping 160 of 160", 158 converged).
Needs gfortran; the tests skip when it is missing.
"""

from __future__ import annotations

import math
import random
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"


@pytest.fixture(scope="module")
def driver(tmp_path_factory) -> Path:
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran not found")
    build = tmp_path_factory.mktemp("ks_select")
    flags = ["-fdefault-real-8", "-fdefault-double-8", "-ffree-line-length-none", "-O0"]
    sources = [
        ROOT / "src" / "argsort.f90",
        ROOT / "src" / "ks_select.f90",
        FORTRAN / "ks_select_driver.f90",
    ]
    for src in sources:
        subprocess.run([gfortran, *flags, "-c", str(src)], cwd=build, check=True)
    exe = build / "driver"
    subprocess.run(
        [gfortran, *flags, "-o", str(exe), *(f"{s.stem}.o" for s in sources)],
        cwd=build,
        check=True,
    )
    return exe


def select(exe: Path, entries, nev: int, delta: float, tol: float = 1.0e-6):
    """entries: list of (re, im, residual). Returns (nsel, n_circle, n_resid, mask)."""
    lines = [f"{len(entries)} {nev} {delta!r} {tol!r}"]
    lines += [f"{re!r} {im!r} {res!r}" for re, im, res in entries]
    run = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", capture_output=True, text=True
    )
    assert run.returncode == 0, run.stdout + run.stderr
    head, mask = run.stdout.split("\n")[:2]
    nsel, n_circle, n_resid = (int(x) for x in head.split())
    return nsel, n_circle, n_resid, [c == "1" for c in mask.strip()]


def pairs_and_reals(n_pairs: int, n_real: int, seed: int):
    """Real Schur order: conjugate pairs in consecutive entries, then real values.

    The two members of a pair carry different residuals, as the Schur vectors do.
    """
    rng = random.Random(seed)
    entries = []
    for _ in range(n_pairs):
        mag = rng.uniform(0.92, 1.0)
        ang = rng.uniform(0.1, 3.0)
        re, im = mag * math.cos(ang), mag * math.sin(ang)
        entries.append((re, im, 10 ** rng.uniform(-8, 0)))
        entries.append((re, -im, 10 ** rng.uniform(-8, 0)))
    for _ in range(n_real):
        entries.append((rng.uniform(0.9, 1.0), 0.0, 10 ** rng.uniform(-8, 0)))
    return entries


def split_pairs(entries, mask):
    """Indices i where entry i and i + 1 are a conjugate pair and only one is selected."""
    return [
        i
        for i in range(len(entries) - 1)
        if entries[i][1] != 0.0
        and entries[i][0] == entries[i + 1][0]
        and entries[i][1] == -entries[i + 1][1]
        and mask[i] != mask[i + 1]
    ]


@pytest.mark.parametrize("nev", [1, 2, 5, 10])
@pytest.mark.parametrize("seed", range(5))
def test_cap_keeps_room_and_never_splits_a_pair(driver, nev, seed):
    # k_dim = 160 with 151 or more values near the unit circle, as in the ellipse report.
    entries = pairs_and_reals(n_pairs=76, n_real=8, seed=seed)
    nsel, _, _, mask = select(driver, entries, nev=nev, delta=0.2)
    assert nsel == sum(mask)
    assert nsel <= len(entries) - nev, "no room left for the Arnoldi steps"
    assert not split_pairs(entries, mask)


def test_all_pairs_cap_is_an_even_count_below_the_limit(driver):
    # 80 pairs, nev = 2: 79 pairs fit under 158. The old code returned 160.
    entries = pairs_and_reals(n_pairs=80, n_real=0, seed=3)
    nsel, n_circle, _, mask = select(driver, entries, nev=2, delta=0.2)
    assert n_circle == 160
    assert nsel == 158
    assert not split_pairs(entries, mask)


def test_the_pair_with_the_largest_residual_is_dropped(driver):
    entries = pairs_and_reals(n_pairs=80, n_real=0, seed=3)
    keys = [max(entries[2 * p][2], entries[2 * p + 1][2]) for p in range(80)]
    worst = keys.index(max(keys))
    _, _, _, mask = select(driver, entries, nev=2, delta=0.2)
    assert not mask[2 * worst] and not mask[2 * worst + 1]
    assert sum(mask) == 158


def test_real_values_keep_the_smallest_residuals(driver):
    residuals = [0.9, 0.1, 0.5, 0.2, 0.8, 0.3, 0.7, 0.4, 0.6, 0.05]
    entries = [(0.95, 0.0, r) for r in residuals]
    nsel, _, _, mask = select(driver, entries, nev=2, delta=0.2)
    assert nsel == 8
    dropped = sorted(range(10), key=lambda i: residuals[i])[8:]
    assert [i for i, m in enumerate(mask) if not m] == sorted(dropped)


def test_selection_below_the_limit_is_unchanged(driver):
    # Two values near the unit circle, one partly converged, the rest far inside.
    entries = [(0.99, 0.0, 0.5), (0.98, 0.0, 0.5), (0.3, 0.0, 1.0e-5)] + [
        (0.1, 0.0, 1.0) for _ in range(7)
    ]
    nsel, n_circle, n_resid, mask = select(driver, entries, nev=2, delta=0.1)
    assert (n_circle, n_resid) == (2, 1)
    # nev + 2 = 4 by magnitude: the three selected plus the largest remaining
    assert nsel == 4
    assert mask[:3] == [True, True, True]


def test_a_conjugate_pair_is_selected_with_its_partner(driver):
    entries = [(0.99, 0.2, 1.0e-9), (0.99, -0.2, 1.0), (0.1, 0.0, 1.0)] + [
        (0.05, 0.0, 1.0) for _ in range(9)
    ]
    _, _, _, mask = select(driver, entries, nev=2, delta=0.05)
    assert mask[0] and mask[1]
