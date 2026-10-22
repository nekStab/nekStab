"""Tests for scripts/compiler_battery.py: the spread, the report and the stale-output rule.

The battery builds and runs stages, which a unit test must not do. These tests cover the
parts that decide what the report says: the spread between profiles, the verdict columns,
the clean environment of a build, and the rule that a source file older than the run is an
error and not a value.
"""

from __future__ import annotations

import importlib
import json
import os
import sys
import time
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"


@pytest.fixture
def battery(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """Import scripts/compiler_battery.py with the scripts folder on the path."""
    monkeypatch.syspath_prepend(str(SCRIPTS))
    return importlib.import_module("compiler_battery")


def make_case(tmp_path: Path, value: float = 0.125) -> Path:
    """Write a stage folder with one reference quantity, residu.dat, a log and DONE."""
    case = tmp_path / "case"
    (case / "ref").mkdir(parents=True)
    (case / "residu.dat").write_text(f"1 0.5\n2 {value}\n")
    (case / "logfile").write_text("total elapsed time : 1.0E+00 sec\n")
    (case / "DONE").write_text("0\n")
    quantity = {"source": "residu.dat:col2:last", "value": 0.125, "tol_rel": 1e-3}
    (case / "ref" / "reference.json").write_text(
        json.dumps({"case": "test", "quantities": {"last_residual": quantity}})
    )
    return case


def test_spread_of_equal_values_is_zero(battery: ModuleType) -> None:
    assert battery.spread([2.0, 2.0, 2.0]) == (0.0, 0.0)


def test_spread_is_relative_to_the_largest_magnitude(battery: ModuleType) -> None:
    absolute, relative = battery.spread([-1.0, 3.0])
    assert absolute == 4.0
    assert relative == pytest.approx(4.0 / 3.0)


def test_spread_of_zeros_does_not_divide_by_zero(battery: ModuleType) -> None:
    assert battery.spread([0.0, 0.0]) == (0.0, 0.0)


def test_slug_has_no_path_separator(battery: ModuleType) -> None:
    assert battery.slug("cylinder_re100/210_baseflow_newton/dyn") == (
        "cylinder_re100-210_baseflow_newton-dyn"
    )


def test_environment_drops_inherited_compiler_variables(
    battery: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("NEKSTAB_FC", "ifx")
    monkeypatch.setenv("CI", "1")
    env = battery.environment(battery.Profile("test", {"NEKSTAB_OPT": "-O0"}))
    assert "NEKSTAB_FC" not in env
    assert "CI" not in env
    assert env["NEKSTAB_OPT"] == "-O0"


def test_every_stage_exists_in_the_examples(battery: ModuleType) -> None:
    for stage in battery.STAGES:
        folder = REPO / "example" / stage.key
        assert (folder / "ref" / "reference.json").is_file(), stage.key
        assert (folder / "run.local.slurm").is_file(), stage.key
        assert (folder / f"{stage.session}.par").is_file(), stage.key


def test_fresh_output_gives_a_value_and_a_passing_check(
    battery: ModuleType, tmp_path: Path
) -> None:
    case = make_case(tmp_path)
    stage = battery.Stage("test/stage", "1cyl", 1, 1)
    result = battery.Result("test", stage.key)
    battery.read_fingerprint(result, case, stage, time.time() - 100.0)
    assert result.errors == {}
    assert result.values == {"last_residual": 0.125}
    assert result.checks == {"last_residual": True}
    assert result.run_gate == ""


def test_output_older_than_the_run_is_an_error(
    battery: ModuleType, tmp_path: Path
) -> None:
    """A tracked copy of an earlier output must not pass as the result of this run."""
    case = make_case(tmp_path)
    old = time.time() - 1000.0
    os.utime(case / "residu.dat", (old, old))
    stage = battery.Stage("test/stage", "1cyl", 1, 1)
    result = battery.Result("test", stage.key)
    battery.read_fingerprint(result, case, stage, time.time() - 100.0)
    assert "last_residual" in result.errors
    assert result.values == {}
    assert result.checks == {}


def test_report_lists_spread_and_verdicts(
    battery: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(battery, "source_state", lambda: "HEAD test")
    stage = battery.Stage("test/stage", "1cyl", 4, 1)
    rows = (
        ("good", 1.0, True, ""),
        ("shifted", 1.5, False, ""),
        ("broken", 0.0, True, "build"),
    )
    for name, value, check, failed in rows:
        folder = tmp_path / name / battery.slug(stage.key)
        folder.mkdir(parents=True)
        result = battery.Result(name, stage.key, failed=failed)
        if not failed:
            result.values["q"] = value
            result.checks["q"] = check
        battery.save(folder, result)
    text = battery.report(tmp_path, [r[0] for r in rows], [stage])
    assert "## test/stage (4 ranks)" in text
    assert "| good |" in text and "| pass |" in text
    assert "| shifted |" in text and "| FAIL |" in text
    assert "build failed" in text
    assert "5.00e-01 | 3.33e-01" in text


def test_main_rejects_an_unknown_profile(
    battery: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(
        sys, "argv", ["compiler_battery.py", str(tmp_path), "--profiles", "nope"]
    )
    assert battery.main() == 2
    assert "unknown profile or stage" in capsys.readouterr().err
