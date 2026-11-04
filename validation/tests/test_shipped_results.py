"""Each field and each figure at a stage root is listed in example/shipped_results.tsv with its md5.

The rules for what a stage ships are in example/REFERENCE_RESULTS.md. A file that a run overwrote fails this test
until `git restore -- <file>` brings the shipped file back, so run the test on a clean checkout.
"""

from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "example" / "shipped_results.tsv"

ROWS = list(csv.DictReader(MANIFEST.open(), delimiter="\t"))


def _tracked() -> list[str]:
    # An unpacked source has no .git: every file under example/ is then a shipped file.
    if not (REPO / ".git").exists():
        return [p.relative_to(REPO).as_posix() for p in (REPO / "example").rglob("*") if p.is_file()]
    return subprocess.run(["git", "-C", str(REPO), "ls-files", "example"], capture_output=True, text=True, check=True).stdout.split("\n")


TRACKED = _tracked()


def test_manifest_has_rows() -> None:
    assert len(ROWS) > 500


def test_each_listed_file_has_its_md5() -> None:
    wrong = [r["path"] for r in ROWS if hashlib.md5((REPO / "example" / r["path"]).read_bytes()).hexdigest() != r["md5"]]
    assert not wrong, wrong


def test_each_field_and_root_figure_is_listed() -> None:
    listed = {"example/" + r["path"] for r in ROWS}
    shipped = [t for t in TRACKED if re.search(r"\.f\d{5}$", t) or (re.search(r"\.(png|gif)$", t) and "/ref/" not in t)]
    missing = [t for t in shipped if t not in listed]
    assert not missing, missing


def _copied_stage(tmp_path: Path) -> Path:
    """A copy of cylinder_re100/310 direct outside the repository: its reference and shipped text files."""
    name = "cylinder_re100/310_stability_direct/direct"
    stage = tmp_path / "copy"
    shutil.copytree(REPO / "example" / name / "ref", stage / "ref")
    for row in ROWS:
        if row["path"].startswith(name + "/") and row["path"].endswith(".dat"):
            shutil.copyfile(REPO / "example" / row["path"], stage / row["path"].rsplit("/", 1)[1])
    return stage


def _check(stage: Path, from_stage: bool) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(REPO / "scripts" / "check_against_ref.py"), "." if from_stage else str(stage)],
        capture_output=True,
        text=True,
        cwd=stage if from_stage else None,
    )


# Every stage README runs the checker from the stage folder with ".", so both forms of the path count.
@pytest.mark.parametrize("from_stage", [False, True])
def test_the_checker_refuses_a_shipped_copy(tmp_path: Path, from_stage: bool) -> None:
    result = _check(_copied_stage(tmp_path), from_stage)
    assert result.returncode != 0
    assert "shipped copy" in result.stdout + result.stderr


def test_the_checker_refuses_a_shipped_copy_in_the_repository() -> None:
    stage = REPO / "example" / "naca0012_re2000" / "310_stability_direct"
    if (stage / "nek5000").exists():
        pytest.fail("the stage holds a build: run this test on a clean checkout")
    assert "shipped copy" in _check(stage, True).stdout


def test_the_checker_accepts_the_same_bytes_written_after_the_build(tmp_path: Path) -> None:
    # A run can repeat a short text file to the last digit.
    stage = _copied_stage(tmp_path)
    exe = stage / "nek5000"
    exe.touch()
    for f in stage.glob("*.dat"):
        os.utime(f, (exe.stat().st_mtime + 10.0,) * 2)
    result = _check(stage, True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_a_shipped_mode_gets_its_coordinates_from_the_stage() -> None:
    # The shipped modes have no coordinates, and every stage README runs plot.py on them.
    sys.path.insert(0, str(REPO / "example"))
    import nekplot

    stage = REPO / "example" / "cylinder_re100" / "310_stability_direct" / "direct"
    x, y, _, _ = nekplot.read_field(stage / "dRe1cyl0.f00001")
    assert x.max() - x.min() > 10.0 and y.max() - y.min() > 10.0


def test_a_field_without_any_coordinates_is_refused(tmp_path: Path) -> None:
    sys.path.insert(0, str(REPO / "example"))
    import nekplot

    source = REPO / "example" / "cylinder_re100" / "310_stability_direct" / "direct" / "dRe1cyl0.f00001"
    alone = tmp_path / source.name
    shutil.copyfile(source, alone)
    with pytest.raises(FileNotFoundError, match="no coordinates"):
        nekplot.read_field(alone)
