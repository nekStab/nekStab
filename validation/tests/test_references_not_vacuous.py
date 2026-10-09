"""Negative control for the example references.

A reference that passes on a stage that never ran checks nothing. The fixture exports the
tracked files of the examples from HEAD into a temporary folder. No stage has run there, so
no output file exists. Every quantity of every ref/reference.json must fail on it. A quantity
that passes reads a tracked input or a tracked copy of an old output.

The control reads HEAD, not the working tree. Commit a change to a stage before it counts.
"""

from __future__ import annotations

import importlib
import json
import subprocess
from pathlib import Path
from types import ModuleType

import pytest

pytest.importorskip("pymech")

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
SUFFIX = "/ref/reference.json"


def tracked_references() -> list[str]:
    """Return the repository paths of the tracked reference files of the examples."""
    if not (REPO / ".git").exists():
        return []
    listing = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z", "--", "example"],
        capture_output=True,
        check=True,
    ).stdout.decode()
    return sorted(p for p in listing.split("\0") if p.endswith(SUFFIX))


REFERENCES = tracked_references()


@pytest.fixture(scope="module")
def unrun_examples(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Export the tracked files of the examples from HEAD. Nothing has run in them."""
    dest = tmp_path_factory.mktemp("unrun")
    archive = subprocess.Popen(
        ["git", "-C", str(REPO), "archive", "HEAD", "example"], stdout=subprocess.PIPE
    )
    subprocess.run(["tar", "-x", "-C", str(dest)], stdin=archive.stdout, check=True)
    if archive.wait() != 0:
        raise RuntimeError("git archive failed")
    return dest


@pytest.fixture
def checker(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """Import scripts/check_against_ref.py."""
    monkeypatch.syspath_prepend(str(SCRIPTS))
    return importlib.import_module("check_against_ref")


def passing_quantities(case: Path, checker: ModuleType) -> list[str]:
    """Return the names of the quantities of a stage folder that pass the reference check."""
    reference = json.loads((case / "ref" / "reference.json").read_text())
    passing = []
    for name, spec in reference.get("quantities", {}).items():
        try:
            value = checker.extract(case, str(spec["source"]))
        except (FileNotFoundError, ValueError):
            continue
        if checker.check_quantity(name, spec, value)[0]:
            passing.append(name)
    return passing


def write_stage(case: Path, with_output: bool) -> None:
    """Write a stage folder with one quantity that reads the last value of residu.dat."""
    (case / "ref").mkdir(parents=True)
    spec = {"last": {"source": "residu.dat:col2:last", "value": 0.125, "tol_rel": 1e-3}}
    (case / "ref" / "reference.json").write_text(json.dumps({"quantities": spec}))
    if with_output:
        (case / "residu.dat").write_text("1 0.5\n2 0.125\n")


def test_control_flags_a_quantity_that_reads_a_present_file(
    tmp_path: Path, checker: ModuleType
) -> None:
    write_stage(tmp_path, with_output=True)
    assert passing_quantities(tmp_path, checker) == ["last"]


def test_control_ignores_a_quantity_whose_output_is_missing(
    tmp_path: Path, checker: ModuleType
) -> None:
    write_stage(tmp_path, with_output=False)
    assert passing_quantities(tmp_path, checker) == []


@pytest.mark.parametrize("reference", REFERENCES)
def test_no_quantity_passes_without_a_run(
    reference: str, unrun_examples: Path, checker: ModuleType
) -> None:
    case = unrun_examples / reference.removesuffix(SUFFIX)
    quantities = json.loads((case / "ref" / "reference.json").read_text()).get(
        "quantities", {}
    )
    assert quantities, f"{reference} lists no quantity, so it checks nothing"
    passing = passing_quantities(case, checker)
    assert not passing, f"{passing} pass on a stage that never ran"


@pytest.mark.parametrize("reference", REFERENCES)
def test_every_quantity_has_a_value_and_a_tolerance(reference: str) -> None:
    quantities = json.loads((REPO / reference).read_text()).get("quantities", {})
    for name, spec in quantities.items():
        assert isinstance(spec.get("value"), int | float), (
            f"{name} has no numeric value"
        )
        assert "source" in spec, f"{name} has no source"
        assert ("tol_rel" in spec) != ("tol_abs" in spec), f"{name} needs one tolerance"
