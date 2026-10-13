"""Ritz value selection of the Krylov-Schur restart (src/ks_select.f90).

The module is compiled alone (it needs no Nek5000) and run once per case by
fortran/ks_select_driver.f90. The cases are built here, not read from a run.

An earlier version removed single values by residual and then completed the
conjugate pairs. The completion returned the count to n, no Arnoldi step
followed, and the residual test of krylov_schur reported n - 2 converged
eigenvalues (an ellipse with k_dim = 160: "keeping 160 of 160", 158 converged).
The export of the eigenmodes (outpost_ks) must visit the converged Ritz values, which
are not always the first entries of the list sorted by modulus (ks_converged_indices).
The run stops when the schur_tgt values at the top of that list have converged
(ks_leading_converged), not when any schur_tgt values have.
A residual below the tolerance does not bound the error of the eigenvalue of a non-normal
operator, so the run also waits until the leading values stand still (ks_leading_settled).
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


@pytest.fixture(scope="module")
def converged_driver(tmp_path_factory) -> Path:
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran not found")
    build = tmp_path_factory.mktemp("ks_converged")
    flags = ["-fdefault-real-8", "-fdefault-double-8", "-ffree-line-length-none", "-O0"]
    sources = [
        ROOT / "src" / "argsort.f90",
        ROOT / "src" / "ks_select.f90",
        FORTRAN / "ks_converged_driver.f90",
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


@pytest.fixture(scope="module")
def settled_driver(tmp_path_factory) -> Path:
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran not found")
    build = tmp_path_factory.mktemp("ks_settled")
    flags = ["-fdefault-real-8", "-fdefault-double-8", "-ffree-line-length-none", "-O0"]
    sources = [
        ROOT / "src" / "argsort.f90",
        ROOT / "src" / "ks_select.f90",
        FORTRAN / "ks_settled_driver.f90",
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


def settled(exe: Path, steps, nev: int = 1, nstill: int = 3, tol: float = 1.0e-6):
    """steps: one list of (re, im, residual) per Arnoldi step. Returns [(settled, residual-only)]."""
    lines = [f"{nev} {nstill} {tol!r}"]
    for entries in steps:
        lines.append(str(len(entries)))
        lines += [f"{re!r} {im!r} {res!r}" for re, im, res in entries]
    run = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", capture_output=True, text=True
    )
    assert run.returncode == 0, run.stdout + run.stderr
    return [tuple(int(x) for x in line.split()) for line in run.stdout.strip().split("\n")]


def converged(exe: Path, entries, tol: float = 1.0e-6):
    """entries: list of (re, im, residual). Returns the 1-based positions of the converged ones."""
    lines = [f"{len(entries)} {tol!r}"]
    lines += [f"{re!r} {im!r} {res!r}" for re, im, res in entries]
    run = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", capture_output=True, text=True
    )
    assert run.returncode == 0, run.stdout + run.stderr
    head, positions = run.stdout.split("\n")[:2]
    assert int(head) == len(positions.split())
    return [int(x) for x in positions.split()]


def leading(exe: Path, entries, tol: float = 1.0e-6) -> int:
    """entries: list of (re, im, residual). Returns the count of ks_leading_converged."""
    lines = [f"{len(entries)} {tol!r}"]
    lines += [f"{re!r} {im!r} {res!r}" for re, im, res in entries]
    run = subprocess.run(
        [str(exe)], input="\n".join(lines) + "\n", capture_output=True, text=True
    )
    assert run.returncode == 0, run.stdout + run.stderr
    return int(run.stdout.split("\n")[2])


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


def test_converged_values_below_an_unconverged_one_are_found(converged_driver):
    # cylinder_re180, Floquet: the phase mode, two unconverged real multipliers of larger
    # modulus than the converged pair 0.233 +- 0.601i (the pair is far below in the list).
    entries = [
        (1.000018, 0.0, 7.7e-12),
        (0.880352, 0.0, 8.8e-3),
        (0.872416, 0.0, 8.3e-3),
        (0.858055, 0.0, 5.2e-4),
        (0.233106, 0.600789, 1.0e-9),
        (0.233106, -0.600789, 1.0e-9),
        (0.1, 0.0, 1.0),
    ]
    assert converged(converged_driver, entries) == [1, 5, 6]


def test_the_first_entries_are_returned_when_they_are_the_converged_ones(converged_driver):
    entries = [(0.99, 0.2, 1.0e-9), (0.99, -0.2, 1.0e-9), (0.8, 0.0, 1.0e-3), (0.5, 0.0, 1.0)]
    assert converged(converged_driver, entries) == [1, 2]


def test_no_value_converged(converged_driver):
    entries = [(0.99, 0.0, 1.0e-2), (0.5, 0.0, 1.0)]
    assert converged(converged_driver, entries) == []


def test_the_tolerance_is_relative_above_modulus_one_and_absolute_below(converged_driver):
    # |lambda| = 3: the bound is tol*3. |lambda| = 0.5: the bound is tol*1, not tol*0.5.
    entries = [(3.0, 0.0, 2.9e-6), (3.0, 0.0, 3.1e-6), (0.5, 0.0, 0.9e-6), (0.5, 0.0, 1.1e-6)]
    assert converged(converged_driver, entries) == [1, 3]


def test_the_leading_count_stops_at_the_first_unconverged_value(converged_driver):
    # A user's cylinder at Re = 42 (Spectre_Hd.dat, multipliers for T = 1): two slow real
    # values lead the list unconverged, the shedding pair below them has converged. A count
    # over the whole list reached schur_tgt = 2 and the run stopped there.
    entries = [
        (0.9863590, 0.0, 6.106864e-3),
        (0.9852301, 0.0, 1.192071e-2),
        (0.7315887, 0.6504731, 2.726931e-9),
        (0.7315887, -0.6504731, 2.726931e-9),
        (0.8977830, 0.2856916, 4.549857e-2),
    ]
    assert converged(converged_driver, entries) == [3, 4]
    assert leading(converged_driver, entries) == 0


def test_a_target_larger_than_the_converged_top_does_not_stop(converged_driver):
    # lid_driven_re3600 310 (Spectre_Hd.dat of the release check): schur_tgt = 4, the real
    # mode and the pair have converged, row 4 has not. 39 values had converged in all.
    entries = [
        (0.9833477, 0.0, 5.303505e-7),
        (0.6872962, -0.6686519, 8.804929e-13),
        (0.6872962, 0.6686519, 8.804929e-13),
        (0.9518588, 0.0, 1.173250e-2),
        (0.9433427, 0.0, 4.502671e-3),
        (0.0088639, -0.9344230, 1.387086e-15),
    ]
    assert leading(converged_driver, entries) == 3


def test_a_converged_pair_at_the_top_counts_as_two(converged_driver):
    entries = [(0.99, 0.2, 1.0e-9), (0.99, -0.2, 1.0e-9), (0.8, 0.0, 1.0e-3), (0.5, 0.0, 1.0e-12)]
    assert leading(converged_driver, entries) == 2


def test_the_leading_count_is_the_whole_list_when_all_converged(converged_driver):
    entries = [(0.99, 0.0, 1.0e-9), (0.9, 0.1, 1.0e-8), (0.9, -0.1, 1.0e-8)]
    assert leading(converged_driver, entries) == 3


def test_a_nan_residual_is_not_converged(converged_driver):
    entries = [(0.99, 0.0, 1.0e-9), (0.9, 0.0, float("nan")), (0.8, 0.0, 1.0e-9)]
    assert leading(converged_driver, entries) == 1
    assert converged(converged_driver, entries) == [1, 3]


def test_a_value_that_stands_still_settles_after_nstill_moves(settled_driver):
    # Residual below the tolerance from step 1. No history at step 1, then three steps of rest.
    steps = [[(1.16, 0.0, 1.0e-9), (0.5, 0.0, 1.0)] for _ in range(6)]
    counts = settled(settled_driver, steps, nev=1, nstill=3)
    assert [c[0] for c in counts] == [0, 0, 0, 1, 1, 1]
    assert [c[1] for c in counts] == [1] * 6


def test_a_value_that_still_moves_does_not_settle_although_its_residual_passes(settled_driver):
    # The forced jet: residual 1.1e-6 (passes 1e-6 * 1.16) while the value moves by 2e-3 per step.
    # The residual test alone says 1 from the first step; the settled count stays 0.
    steps = [[(1.16 + 2.0e-3 * 0.7**k, 0.0, 1.0e-7)] for k in range(8)]
    counts = settled(settled_driver, steps, nev=1, nstill=3)
    assert all(c == (0, 1) for c in counts)


def test_a_value_settles_once_its_steps_are_below_the_tolerance(settled_driver):
    # Moves of 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9: below 1.16e-6 from the move of 1e-6 on.
    # The moves 1e-6 (not below 1.16e-6? it is: 1e-6 < 1.16e-6), 1e-7, 1e-8 are three steps of rest.
    values = [1.16]
    for k in range(3, 10):
        values.append(values[-1] + 10.0 ** (-k))
    steps = [[(v, 0.0, 1.0e-9)] for v in values]
    counts = settled(settled_driver, steps, nev=1, nstill=3)
    # step 1 has no history; steps 2 and 3 move by 1e-3 and 1e-4; the move 1e-6 enters at step 5
    assert [c[0] for c in counts] == [0, 0, 0, 0, 0, 0, 1, 1]


def test_a_jump_resets_the_count(settled_driver):
    values = [1.16] * 5 + [1.17] + [1.17] * 4
    steps = [[(v, 0.0, 1.0e-9)] for v in values]
    counts = settled(settled_driver, steps, nev=1, nstill=3)
    assert [c[0] for c in counts] == [0, 0, 0, 1, 1, 0, 0, 0, 1, 1]


def test_a_conjugate_pair_that_changes_places_is_still_stationary(settled_driver):
    # The members of a pair may come in either order from one step to the next.
    first = [(0.9, 0.3, 1.0e-9), (0.9, -0.3, 1.0e-9)]
    second = [(0.9, -0.3, 1.0e-9), (0.9, 0.3, 1.0e-9)]
    steps = [first, second, first, second, first]
    counts = settled(settled_driver, steps, nev=2, nstill=3)
    assert [c[0] for c in counts] == [0, 0, 0, 2, 2]


def test_an_unconverged_value_does_not_settle_however_still_it_is(settled_driver):
    steps = [[(1.16, 0.0, 1.0e-3)] for _ in range(6)]
    counts = settled(settled_driver, steps, nev=1, nstill=3)
    assert all(c == (0, 0) for c in counts)


def test_the_target_needs_every_leading_entry_to_stand_still(settled_driver):
    # nev = 2: the second entry keeps moving, the first does not.
    steps = [[(1.16, 0.0, 1.0e-9), (1.0 - 0.1 * 0.5**k, 0.0, 1.0e-9)] for k in range(7)]
    counts = settled(settled_driver, steps, nev=2, nstill=3)
    assert all(c[0] <= 1 for c in counts)
    assert counts[-1][0] <= 1


def test_fewer_values_than_targets_cannot_reach_the_target(settled_driver):
    # One Ritz value at step 1 to 6 and nev = 3: the count can reach 1 but never 3.
    steps = [[(1.16, 0.0, 1.0e-9)] for _ in range(6)]
    counts = settled(settled_driver, steps, nev=3, nstill=3)
    assert [c[0] for c in counts] == [0, 0, 0, 1, 1, 1]
    assert all(c[0] < 3 for c in counts)


def test_a_nan_value_is_not_stationary(settled_driver):
    steps = [[(1.16, 0.0, 1.0e-9)], [(float("nan"), 0.0, 1.0e-9)], [(1.16, 0.0, 1.0e-9)]]
    counts = settled(settled_driver, steps, nev=1, nstill=1)
    assert [c[0] for c in counts] == [0, 0, 0]
