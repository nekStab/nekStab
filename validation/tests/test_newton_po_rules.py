"""Check the Newton periodic-orbit step, CFL warning and filter warning.

The arithmetic is in src/newton_po_rules.f90. A kept step count gives
dt = T/N. The unfrozen branch also ends at dt = T/N, with N from the
provisional step. The CFL warning is strict: 1.5 times the target does
not warn. An explicit filter with a positive weight warns only for a
periodic-orbit Newton mode. Needs gfortran. The tests skip when it is missing.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"

# filterType: 0 none, 1 explicit, 2 relaxation. Weight 0 and a positive weight.
FILTER_TYPES = (0, 1, 2)
FILTER_WEIGHTS = (0.0, 0.01)
# (is_po, is_po_t). Both orbit modes warn. A steady mode does not.
ORBIT_MODES = ((1, 0), (0, 1), (0, 0))


@pytest.fixture(scope="module")
def driver(build_driver) -> Path:
    return build_driver(
        "newton_po_rules",
        [
            ROOT / "src" / "newton_po_rules.f90",
            FORTRAN / "newton_po_rules_driver.f90",
        ],
    )


def ask(exe: Path, text: str) -> list[str]:
    run = subprocess.run([str(exe)], input=text, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    return run.stdout.split()


def test_kept_count_sets_dt_to_period_over_n(driver):
    # dt_in = 0.3 would give another count. The kept count wins.
    dt, nsteps = ask(driver, "step\n8 4 0.3\n")
    assert int(nsteps) == 4
    assert float(dt) == 2.0

    dt, nsteps = ask(driver, "step\n1 8 0.3\n")
    assert int(nsteps) == 8
    assert float(dt) == 0.125


def test_unfrozen_branch_sets_dt_to_period_over_count(driver):
    # The provisional step 0.25 already divides the period.
    dt, nsteps = ask(driver, "step\n1 0 0.25\n")
    assert int(nsteps) == 4
    assert float(dt) == 0.25

    # ceiling(1/0.3) = 4, so the step becomes T/N and not 0.3.
    dt, nsteps = ask(driver, "step\n1 0 0.3\n")
    assert int(nsteps) == 4
    assert float(dt) == 0.25


@pytest.mark.parametrize(
    ("n_kept", "target", "cfl", "warns"),
    [
        (4, 2.0, 3.0, False),  # exactly 1.5 x target
        (4, 2.0, 3.2, True),  # above 1.5 and below 2.0
        (4, 0.0, 10.0, False),  # target 0
        (0, 2.0, 10.0, False),  # no kept count
    ],
)
def test_cfl_warning_predicate(driver, n_kept, target, cfl, warns):
    (flag,) = ask(driver, f"cfl\n{n_kept} {target} {cfl}\n")
    assert int(flag) == int(warns)


@pytest.mark.parametrize("is_po,is_po_t", ORBIT_MODES)
@pytest.mark.parametrize("filter_type", FILTER_TYPES)
@pytest.mark.parametrize("weight", FILTER_WEIGHTS)
def test_explicit_filter_warning_predicate(driver, is_po, is_po_t, filter_type, weight):
    orbit = is_po == 1 or is_po_t == 1
    expects = orbit and filter_type == 1 and weight > 0.0
    (flag,) = ask(driver, f"filter\n{is_po} {is_po_t} {filter_type} {weight}\n")
    assert int(flag) == int(expects)


@pytest.mark.parametrize(
    ("iteration", "is_po", "is_po_t", "nsteps", "stored", "kept"),
    [
        (1, 1, 0, 1212, 0, 1212),
        (1, 0, 1, 800, 0, 800),
        (1, 0, 0, 100, 0, 0),
        (2, 1, 0, 1300, 1212, 1212),
        (5, 0, 1, 900, 800, 800),
        (1, 1, 1, 50, 7, 50),
    ],
)
def test_kept_count_is_set_at_iteration_1_for_an_orbit(
    driver, iteration, is_po, is_po_t, nsteps, stored, kept
):
    (out,) = ask(
        driver, f"keep\n{iteration} {is_po} {is_po_t} {nsteps} {stored}\n"
    )
    assert int(out) == kept


def test_kept_count_is_cleared_after_the_stage(driver):
    (out,) = ask(driver, "clear\n")
    assert int(out) == 0
