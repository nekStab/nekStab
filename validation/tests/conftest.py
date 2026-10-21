"""Shared fixture: build a gfortran test driver from Fortran sources.

The drivers compile a nekStab module alone (no Nek5000) with a small main program.
The tests that use the fixture skip when gfortran is missing.
"""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path

import pytest

FLAGS = ["-fdefault-real-8", "-fdefault-double-8", "-ffree-line-length-none", "-O0"]


@pytest.fixture(scope="session")
def build_driver(
    tmp_path_factory: pytest.TempPathFactory,
) -> Callable[[str, Sequence[Path]], Path]:
    """Return a function that compiles the sources in order and links one executable."""
    gfortran = shutil.which("gfortran")

    def build(name: str, sources: Sequence[Path]) -> Path:
        if gfortran is None:
            pytest.skip("gfortran not found")
        # The sources compile in the given order: each module comes before its users.
        work = tmp_path_factory.mktemp(name)
        for src in sources:
            subprocess.run([gfortran, *FLAGS, "-c", str(src)], cwd=work, check=True)
        exe = work / "driver"
        subprocess.run(
            [gfortran, *FLAGS, "-o", str(exe), *(f"{s.stem}.o" for s in sources)],
            cwd=work,
            check=True,
        )
        return exe

    return build
