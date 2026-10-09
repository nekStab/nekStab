#!/usr/bin/env python3
"""Run example stages under several compiler profiles and compare the results.

For each pair of profile and stage the script exports the tracked files of the
stage from HEAD, builds the case with bin/makeneks, runs the run script of the
stage, and stores a fingerprint. The fingerprint holds every quantity of
ref/reference.json at the precision of the output files, and a few extra
sources of the stage. The report lists the value of each quantity under each
profile, the spread between the profiles, and the result of the reference check.

Usage:
    python3 scripts/compiler_battery.py WORKDIR [--profiles NAME,...] [--stages KEY,...]
    python3 scripts/compiler_battery.py WORKDIR --report-only

A folder WORKDIR/<profile>/<stage> that holds a result.json is not run again.
Delete the folder to repeat the pair. Run the script under jobscope, because it
starts builds and MPI jobs. The build uses the source files of the working tree
and the stage files of HEAD; the report states if src or bin differs from HEAD.
"""

import argparse
import json
import math
import os
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import check_against_ref as car

REPO = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Profile:
    """Compiler selection and flags of one build."""

    name: str
    env: dict[str, str]


@dataclass(frozen=True)
class Stage:
    """One example stage of the battery.

    `needs` lists stages that are exported next to it. `copies` lists glob
    patterns, relative to the stage folder, of files copied into the folder
    before the run (the inputs that the README of the stage tells the user to copy).
    `extra` maps names to sources that are added to the fingerprint.
    """

    key: str
    session: str
    ranks: int
    minutes: int
    needs: tuple[str, ...] = ()
    copies: tuple[str, ...] = ()
    extra: tuple[tuple[str, str], ...] = ()


@dataclass
class Result:
    """Outcome of one profile and stage, stored as result.json."""

    profile: str
    stage: str
    build_s: float = 0.0
    run_s: float = 0.0
    run_rc: int | None = None
    failed: str = ""
    run_gate: str = ""
    values: dict[str, float] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)
    checks: dict[str, bool] = field(default_factory=dict)


TRAPS_GNU = "-g -ffpe-trap=invalid,zero,overflow"
SNAN_GNU = "-g -finit-real=snan"
# The battery has no Intel profile with traps or with signaling NaNs. With -fpe0 a
# run stops in MKL (mkl_lapack_dlaq5, called by dsygv in the Nek5000 setup of the
# Schwarz smoother) with a division by zero. The Intel flag -init=snan,arrays also
# reaches static arrays, and Nek5000 relies on zeros there.

PROFILES = (
    Profile("gnu-default", {"NEKSTAB_FC": "gcc"}),
    Profile("gnu-portable", {"NEKSTAB_FC": "gcc", "CI": "1"}),
    Profile("gnu-O3", {"NEKSTAB_FC": "gcc", "NEKSTAB_OPT": "-O3"}),
    Profile("gnu-O0", {"NEKSTAB_FC": "gcc", "NEKSTAB_OPT": "-O0"}),
    Profile(
        "gnu-traps",
        {"NEKSTAB_FC": "gcc", "NEKSTAB_OPT": "-O0", "NEKSTAB_EXTRA_FFLAGS": TRAPS_GNU},
    ),
    Profile(
        "gnu-snan",
        {"NEKSTAB_FC": "gcc", "NEKSTAB_OPT": "-O0", "NEKSTAB_EXTRA_FFLAGS": SNAN_GNU},
    ),
    Profile("ifx-default", {"NEKSTAB_FC": "ifx"}),
    Profile("ifx-portable", {"NEKSTAB_FC": "ifx", "CI": "1"}),
    Profile("ifx-O0", {"NEKSTAB_FC": "ifx", "NEKSTAB_OPT": "-O0"}),
)

FLIP_FLOP = "flip_flop_re62"
STAGES = (
    Stage(
        "cylinder_re100/000_dns",
        "1cyl",
        8,
        30,
        extra=(
            ("snapshot_vx_norm", "1cyl0.f00020:field:vx.norm"),
            ("snapshot_vy_norm", "1cyl0.f00020:field:vy.norm"),
            ("snapshot_vy_max", "1cyl0.f00020:field:vy.max"),
        ),
    ),
    Stage(
        "thermosyphon_ra500/210_baseflow_newton",
        "tsyphon",
        4,
        20,
        extra=(
            ("bf_vx_norm", "BF_tsyphon0.f00001:field:vx.norm"),
            ("bf_vy_norm", "BF_tsyphon0.f00001:field:vy.norm"),
            ("bf_t_norm", "BF_tsyphon0.f00001:field:t.norm"),
            ("bf_t_max", "BF_tsyphon0.f00001:field:t.max"),
        ),
    ),
    Stage(
        "thermosyphon_ra500/310_stability_direct/direct",
        "tsyphon",
        4,
        20,
        extra=(
            ("sigma_row2", "Spectre_NSd.dat:col1:row2"),
            ("sigma_row3", "Spectre_NSd.dat:col1:row3"),
            ("omega_row2", "Spectre_NSd.dat:col2:row2"),
        ),
    ),
    Stage("cylinder_re100/120_baseflow_boostconv", "1cyl", 8, 30),
    Stage("cylinder_re100/210_baseflow_newton/dyn", "1cyl", 8, 60),
    Stage(
        f"{FLIP_FLOP}/411_postproc_wavemaker",
        "2cyl",
        8,
        20,
        needs=(
            f"{FLIP_FLOP}/311_stability_direct_floquet",
            f"{FLIP_FLOP}/321_stability_adjoint_floquet",
        ),
        copies=(
            "../311_stability_direct_floquet/dRe2cyl0.f0000[12]",
            "../311_stability_direct_floquet/dIm2cyl0.f0000[12]",
            "../321_stability_adjoint_floquet/aRe2cyl0.f0000[12]",
            "../321_stability_adjoint_floquet/aIm2cyl0.f0000[12]",
        ),
    ),
    Stage(
        f"{FLIP_FLOP}/311_stability_direct_floquet",
        "2cyl",
        16,
        90,
        extra=(
            ("sigma_row3", "Spectre_NSd.dat:col1:row3"),
            ("sigma_row5", "Spectre_NSd.dat:col1:row5"),
            ("omega_row6", "Spectre_NSd.dat:col2:row6"),
        ),
    ),
)

# The Newton solve takes 11 min on 8 ranks and the Floquet eigenvalue solve 20 min
# on 16 ranks. The default selection leaves them out.
LONG_STAGES = (
    "cylinder_re100/210_baseflow_newton/dyn",
    f"{FLIP_FLOP}/311_stability_direct_floquet",
)
DEFAULT_STAGES = tuple(s.key for s in STAGES if s.key not in LONG_STAGES)


def slug(text: str) -> str:
    """Return a folder name for a stage key."""
    return text.replace("/", "-")


def export_stage(stage: Stage, dest: Path) -> Path:
    """Write the tracked files of the stage and its helper stages into dest."""
    paths = [f"example/{key}" for key in (stage.key, *stage.needs)]
    dest.mkdir(parents=True)
    archive = subprocess.Popen(
        ["git", "-C", str(REPO), "archive", "HEAD", *paths], stdout=subprocess.PIPE
    )
    subprocess.run(["tar", "-x", "-C", str(dest)], stdin=archive.stdout, check=True)
    if archive.wait() != 0:
        raise RuntimeError(f"git archive failed for {stage.key}")
    return dest / "example" / stage.key


def environment(profile: Profile) -> dict[str, str]:
    """Return the environment of a build and a run under the profile."""
    env = dict(os.environ)
    for name in ("NEKSTAB_FC", "NEKSTAB_OPT", "NEKSTAB_EXTRA_FFLAGS", "CI"):
        env.pop(name, None)
    env["NEKSTAB_SOURCE_ROOT"] = str(REPO)
    env.setdefault("NEK_SOURCE_ROOT", str(REPO / "Nek5000"))
    env["PATH"] = f"{REPO / 'bin'}:{env['NEK_SOURCE_ROOT']}/bin:{env['PATH']}"
    env.update(profile.env)
    return env


def run_logged(
    command: list[str], cwd: Path, env: dict[str, str], log: Path, minutes: int
) -> int:
    """Run a command with its output in a log file.

    Return the exit code, or -1 after the time limit.
    """
    with log.open("w") as handle:
        proc = subprocess.Popen(
            command,
            cwd=cwd,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            return proc.wait(timeout=60 * minutes)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            return -1


def read_fingerprint(result: Result, case: Path, stage: Stage, started: float) -> None:
    """Store the quantities of the reference and the extra sources of a finished run.

    A source whose file is older than `started` counts as an error: the run did
    not write it, and the file is a tracked copy of an earlier output.
    """
    spec: dict[str, dict[str, object]] = json.loads(
        (case / "ref" / "reference.json").read_text()
    )["quantities"]
    sources = {name: str(q["source"]) for name, q in spec.items()}
    sources.update(dict(stage.extra))
    for name, source in sources.items():
        path = case / source.split(":")[0]
        try:
            if path.exists() and path.stat().st_mtime < started:
                raise ValueError(f"{path.name} is older than the run")
            result.values[name] = car.extract(case, source)
        except (FileNotFoundError, ValueError) as exc:
            result.errors[name] = str(exc)
    for name, q in spec.items():
        if name in result.values:
            result.checks[name] = car.check_quantity(name, q, result.values[name])[0]
    result.run_gate = car.run_failed(case) or ""


def save(folder: Path, result: Result) -> None:
    """Write result.json into the folder of the pair."""
    (folder / "result.json").write_text(json.dumps(asdict(result), indent=1))


def run_pair(profile: Profile, stage: Stage, workdir: Path) -> None:
    """Build and run one stage under one profile and write its result.json."""
    folder = workdir / profile.name / slug(stage.key)
    if (folder / "result.json").exists():
        print(f"skip {profile.name} {stage.key}: result.json exists", flush=True)
        return
    if folder.exists():
        shutil.rmtree(folder)
    case = export_stage(stage, folder)
    env = environment(profile)
    for pattern in stage.copies:
        for source in sorted(case.glob(pattern)):
            shutil.copy(source, case)
    (case / "SESSION.NAME").write_text(f"{stage.session}\n{case}/\n")

    result = Result(profile.name, stage.key)
    t0 = time.monotonic()
    build_rc = run_logged(
        [str(REPO / "bin" / "makeneks"), stage.session],
        case,
        env,
        folder / "build.log",
        30,
    )
    result.build_s = round(time.monotonic() - t0, 1)
    if build_rc != 0 or not (case / "nek5000").exists():
        result.failed = "build"
        save(folder, result)
        print(
            f"{profile.name} {stage.key}: BUILD FAILED, see {folder / 'build.log'}",
            flush=True,
        )
        return

    started = time.time()
    run_env = {**env, "SLURM_NTASKS": str(stage.ranks), "SLURM_JOB_ID": "battery"}
    t0 = time.monotonic()
    result.run_rc = run_logged(
        ["bash", "run.local.slurm"], case, run_env, case / "logfile", stage.minutes
    )
    result.run_s = round(time.monotonic() - t0, 1)
    (case / "DONE").write_text(f"{result.run_rc}\n")
    read_fingerprint(result, case, stage, started)
    save(folder, result)
    clean = result.run_rc == 0 and not result.errors and all(result.checks.values())
    print(
        f"{profile.name} {stage.key}: {'ok' if clean else 'NOT CLEAN'}"
        f" build {result.build_s} s run {result.run_s} s",
        flush=True,
    )


def spread(values: list[float]) -> tuple[float, float]:
    """Return the absolute spread and the spread relative to the largest magnitude."""
    absolute = max(values) - min(values)
    scale = max(abs(v) for v in values)
    return absolute, (absolute / scale if scale > 0.0 else 0.0)


def load_results(workdir: Path, profiles: list[str], stage: Stage) -> list[Result]:
    """Return the result of every profile that has run the stage."""
    found = []
    for name in profiles:
        path = workdir / name / slug(stage.key) / "result.json"
        if path.exists():
            found.append(Result(**json.loads(path.read_text())))
    return found


def source_state() -> str:
    """Return the HEAD commit and the state of src and bin against it."""
    head = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "-C", str(REPO), "status", "--porcelain", "--", "src", "bin"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    return f"HEAD {head}, " + (
        "src and bin differ from it" if dirty else "src and bin as in it"
    )


def stage_tables(stage: Stage, results: list[Result]) -> list[str]:
    """Return the Markdown tables of one stage: run summary and quantities."""
    lines = [f"## {stage.key} ({stage.ranks} ranks)", ""]
    lines += [
        "| profile | build s | run s | exit | run gate | reference check |",
        "|---|---:|---:|---:|---|---|",
    ]
    for res in results:
        if res.failed:
            lines.append(f"| {res.profile} | {res.build_s} | | | build failed | |")
            continue
        verdict = "pass" if res.checks and all(res.checks.values()) else "FAIL"
        lines.append(
            f"| {res.profile} | {res.build_s} | {res.run_s} | {res.run_rc} "
            f"| {res.run_gate or 'clean'} | {verdict} |"
        )
    ran = [r for r in results if not r.failed]
    quantities = sorted({q for r in ran for q in r.values})
    if not quantities:
        return [*lines, ""]
    lines += [
        "",
        "| quantity | " + " | ".join(r.profile for r in ran) + " | spread | relative |",
        "|---|" + "---:|" * (len(ran) + 2),
    ]
    for quantity in quantities:
        row = [r.values.get(quantity) for r in ran]
        have = [v for v in row if v is not None and math.isfinite(v)]
        cells = [f"{v:.10g}" if v is not None else "missing" for v in row]
        tail = " | |"
        if len(have) > 1:
            absolute, relative = spread(have)
            tail = f" | {absolute:.2e} | {relative:.2e} |"
        lines.append(f"| {quantity} | " + " | ".join(cells) + tail)
    return [*lines, ""]


def report(workdir: Path, profiles: list[str], stages: list[Stage]) -> str:
    """Return a Markdown report of every result in the working folder."""
    lines = ["# Compiler battery", "", source_state(), ""]
    for stage in stages:
        results = load_results(workdir, profiles, stage)
        if results:
            lines += stage_tables(stage, results)
    return "\n".join(lines)


def main() -> int:
    """Run the selected pairs and print the report."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--profiles", default=",".join(p.name for p in PROFILES))
    parser.add_argument("--stages", default=",".join(DEFAULT_STAGES))
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()

    by_profile = {p.name: p for p in PROFILES}
    by_stage = {s.key: s for s in STAGES}
    profile_names = args.profiles.split(",")
    stage_keys = args.stages.split(",")
    unknown = [n for n in profile_names if n not in by_profile]
    unknown += [k for k in stage_keys if k not in by_stage]
    if unknown:
        print(f"unknown profile or stage: {unknown}", file=sys.stderr)
        return 2

    args.workdir.mkdir(parents=True, exist_ok=True)
    if not args.report_only:
        for name in profile_names:
            for key in stage_keys:
                run_pair(by_profile[name], by_stage[key], args.workdir)
    text = report(args.workdir, profile_names, [by_stage[k] for k in stage_keys])
    (args.workdir / "report.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
