"""Each field and each figure at a stage root is listed in validation/shipped_results.tsv with its md5.

The rules for what a stage ships are in example/REFERENCE_RESULTS.md. A file that a run overwrote fails this test
until `git checkout -- <stage>` brings the shipped file back, so run the test on a clean checkout.
"""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "validation" / "shipped_results.tsv"

if not (REPO / ".git").exists():
    pytest.skip("needs a git checkout", allow_module_level=True)

ROWS = list(csv.DictReader(MANIFEST.open(), delimiter="\t"))
TRACKED = subprocess.run(["git", "-C", str(REPO), "ls-files", "example"], capture_output=True, text=True, check=True).stdout.split("\n")


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
