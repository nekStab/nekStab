"""Every example stage that uses history points ships its probe list.

Nek5000 hpts reads the probe list from <case>.his and stops with "Cannot open history file in
subroutine hpts()" when the file is missing. *.his is gitignored, so each stage tracks the list as
<case>.probes and bin/mks copies it to <case>.his when the build starts. These tests fail when

- a stage that calls hpts (in its .usr, or by a run script that writes <case>.his) has neither a
  tracked <case>.probes nor a tracked <case>.his;
- the first line of a .probes file is not the number of the coordinate lines after it;
- a .probes file differs from the text that the run script of its stage writes to <case>.his;
- the mks step does not copy the list, overwrites an existing .his, or does not warn when both files
  are missing.

The mks step is a bash function of bin/makeneks (mks is a link to it). The tests cut it out of the
script and run it in a temporary folder, so they need no compiler.
"""

from __future__ import annotations

import codecs
import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
EXAMPLE = REPO / "example"
MAKENEKS = REPO / "bin" / "makeneks"

# A stage with no probe line in its own run script: the probe list is the one of its sibling stage,
# which has the same .usr apart from the mode.
SIBLING = {
    "example/cubic_cavity_re1950/411_postproc_energy_budget": "example/cubic_cavity_re1950/311_stability_direct_floquet",
}

WRITE_LINE = re.compile(r"""^\s*printf\s+(['"])(.*?)\1\s*>\s*(\S+?)\.his\s*$""")


def tracked_files() -> list[str]:
    """Return the files that git tracks and that exist below example/, or all files outside a git checkout."""
    if (REPO / ".git").exists():
        out = subprocess.run(
            ["git", "-C", str(REPO), "ls-files", "example"], capture_output=True, text=True, check=True
        )
        # A file removed from the working tree is still in the index until git rm.
        return [f for f in out.stdout.splitlines() if (REPO / f).is_file()]
    return [p.relative_to(REPO).as_posix() for p in EXAMPLE.rglob("*") if p.is_file()]


FILES = tracked_files()


def folder_files(folder: str) -> list[str]:
    """Return the tracked files that sit directly in a folder."""
    prefix = folder + "/"
    return [f for f in FILES if f.startswith(prefix) and "/" not in f[len(prefix) :]]


def calls_hpts(usr: Path) -> bool:
    """True when the fixed-form source has a call of hpts that is not a comment."""
    for line in usr.read_text(errors="replace").splitlines():
        if line[:1] in "cC*!":
            continue
        if re.search(r"\bcall\s+hpts\b", line.split("!")[0], re.I):
            return True
    return False


def script_probe_text(folder: str) -> str | None:
    """Return the text that the run scripts of a folder write to <case>.his, or None."""
    texts = set()
    for f in folder_files(folder):
        if f.endswith((".slurm", ".sh")):
            for line in (REPO / f).read_text().splitlines():
                m = WRITE_LINE.match(line)
                if m:
                    texts.add(codecs.decode(m.group(2), "unicode_escape"))
    assert len(texts) <= 1, f"{folder}: run scripts write different probe lists: {texts}"
    return texts.pop() if texts else None


def hpts_stages() -> list[str]:
    """Return the stage folders that need a probe list."""
    stages = set()
    for f in FILES:
        if f.endswith(".usr") and calls_hpts(REPO / f):
            stages.add(f.rsplit("/", 1)[0])
    for f in FILES:
        # A stage can reach hpts through a routine of src/ (mode 4.52) and have no call in its .usr.
        if f.endswith((".slurm", ".sh")) and script_probe_text(f.rsplit("/", 1)[0]) is not None:
            stages.add(f.rsplit("/", 1)[0])
    return sorted(stages)


STAGES = hpts_stages()


def probe_files(folder: str) -> list[str]:
    return [f for f in folder_files(folder) if f.endswith(".probes")]


PROBE_STAGES = [s for s in STAGES if probe_files(s)]


def test_stages_were_found() -> None:
    assert len(STAGES) >= 40, f"only {len(STAGES)} stages found; the search is broken"


@pytest.mark.parametrize("stage", STAGES)
def test_stage_has_a_probe_list(stage: str) -> None:
    has_probes = bool(probe_files(stage))
    has_his = any(f.endswith(".his") for f in folder_files(stage))
    assert has_probes or has_his, (
        f"{stage} calls hpts but tracks neither <case>.probes nor <case>.his; "
        "mks cannot write the history file, and hpts stops without it"
    )


@pytest.mark.parametrize("stage", STAGES)
def test_probe_file_format(stage: str) -> None:
    for f in probe_files(stage):
        lines = (REPO / f).read_text().splitlines()
        assert lines, f"{f} is empty"
        assert re.fullmatch(r"\s*\d+\s*", lines[0]), f"{f}: first line {lines[0]!r} is not an integer"
        coords = lines[1:]
        assert len(coords) == int(lines[0]), f"{f}: first line says {lines[0].strip()} probes, {len(coords)} lines follow"
        for c in coords:
            values = c.split()
            assert len(values) >= 2, f"{f}: coordinate line {c!r} has fewer than 2 numbers"
            [float(v) for v in values]


@pytest.mark.parametrize("stage", PROBE_STAGES)
def test_probe_file_equals_the_script(stage: str) -> None:
    files = probe_files(stage)
    text = script_probe_text(stage)
    if text is None:
        text = script_probe_text(SIBLING[stage])
        assert text is not None, f"{stage}: no run script writes a probe list, and no sibling is named"
    for f in files:
        assert (REPO / f).read_text() == text, f"{f} differs from the probe list its run script writes"


def test_probe_file_is_named_after_the_case() -> None:
    """mks looks for <case>.probes, with <case> the name of the .usr file."""
    for stage in STAGES:
        usr = [Path(f).stem for f in folder_files(stage) if f.endswith(".usr")]
        for f in probe_files(stage):
            assert Path(f).stem in usr, f"{f} has no .usr of the same name"


# ---------------------------------------------------------------------------------------------
# The mks step
# ---------------------------------------------------------------------------------------------


def run_mks_step(tmp_path: Path, case: str = "c") -> subprocess.CompletedProcess:
    """Run ensure_probe_file of bin/makeneks in tmp_path."""
    text = MAKENEKS.read_text()
    m = re.search(r"^ensure_probe_file \(\) \{\n.*?^\}\n", text, re.S | re.M)
    assert m, "ensure_probe_file is missing in bin/makeneks"
    script = tmp_path / "step.sh"
    script.write_text(m.group(0) + f'ensure_probe_file "{case}"\n')
    return subprocess.run(["bash", str(script)], cwd=tmp_path, capture_output=True, text=True)


def test_mks_is_a_link_to_makeneks() -> None:
    assert (REPO / "bin" / "mks").resolve() == MAKENEKS.resolve()


def test_mks_calls_the_step_before_the_build() -> None:
    text = MAKENEKS.read_text()
    assert 'ensure_probe_file "$target"' in text


def test_mks_copies_the_probe_list(tmp_path: Path) -> None:
    (tmp_path / "c.usr").write_text("      call hpts\n")
    (tmp_path / "c.probes").write_text("2\n1 0 0\n2 0 0\n")
    r = run_mks_step(tmp_path)
    assert r.returncode == 0
    assert (tmp_path / "c.his").read_text() == "2\n1 0 0\n2 0 0\n"
    assert "c.his written from c.probes" in r.stdout
    assert r.stderr == ""


def test_mks_stops_when_the_copy_fails(tmp_path: Path) -> None:
    (tmp_path / "c.usr").write_text("      call hpts\n")
    (tmp_path / "c.probes").write_text("1\n5 0 0\n")
    # c.his is a link to a missing folder: [ -f ] is false and cp fails, also for root
    (tmp_path / "c.his").symlink_to(tmp_path / "missing" / "c.his")
    r = run_mks_step(tmp_path)
    assert r.returncode != 0
    assert "c.his written" not in r.stdout
    assert "cannot copy c.probes" in r.stderr
    assert not (tmp_path / "c.his").exists()


def test_mks_does_not_overwrite_a_history(tmp_path: Path) -> None:
    history = "1\n5 0 0\n  1.0E+00  1.0E-01  2.0E-01  0.0E+00\n"
    (tmp_path / "c.usr").write_text("      call hpts\n")
    (tmp_path / "c.probes").write_text("1\n5 0 0\n")
    (tmp_path / "c.his").write_text(history)
    r = run_mks_step(tmp_path)
    assert r.returncode == 0
    assert (tmp_path / "c.his").read_text() == history
    assert "written" not in r.stdout


def test_mks_warns_when_both_files_are_missing(tmp_path: Path) -> None:
    (tmp_path / "c.usr").write_text("      subroutine userchk\n      call hpts\n")
    r = run_mks_step(tmp_path)
    assert r.returncode == 0
    assert not (tmp_path / "c.his").exists()
    assert "c.probes" in r.stderr and "number of probes" in r.stderr


@pytest.mark.parametrize(
    "usr",
    [
        "c     call hpts\n",
        "      ! call hpts\n",
        "!      call hpts\n",
        "      call hpts_other_name\n",
        "      x = 1 ! call hpts\n",
        "      call nekStab\n",
    ],
)
def test_mks_is_quiet_without_a_call_of_hpts(tmp_path: Path, usr: str) -> None:
    (tmp_path / "c.usr").write_text(usr)
    r = run_mks_step(tmp_path)
    assert r.returncode == 0
    assert r.stdout == "" and r.stderr == ""
    assert not (tmp_path / "c.his").exists()


def test_mks_finds_a_call_after_an_if(tmp_path: Path) -> None:
    (tmp_path / "c.usr").write_text("      if (istep.gt.0) call HPTS\n")
    r = run_mks_step(tmp_path)
    assert "c.probes" in r.stderr


def _slurm_scripts():
    return sorted(p for p in (REPO / "example").rglob("run.local.slurm") if (p.parent / "SIZE").exists())


@pytest.mark.parametrize("script", _slurm_scripts(), ids=lambda p: str(p.parent.relative_to(REPO / "example")))
def test_rank_default_of_the_script_is_not_below_lpmin(script):
    # Without Slurm, SLURM_NTASKS is unset and the script's default decides the
    # rank count. Below lpmin, Nek5000 stops at start-up with a message about the
    # array sizes, which tells a new user nothing about the script.
    lpmin = int(re.search(r"lpmin\s*=\s*(\d+)", (script.parent / "SIZE").read_text()).group(1))
    for default in re.findall(r"SLURM_NTASKS:-(\d+)", script.read_text()):
        assert int(default) >= lpmin
