#!/usr/bin/env python3
from __future__ import annotations

import argparse
import configparser
import glob
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path


SCHEMA_VERSION = "0.1.0"
ARTIFACT_GLOBS = ("BF_*.f*", "*.f0????", "residu*.dat", "logfile", "*.png", "*.mp4")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def utc_mtime(path: Path) -> str:
    return (
        datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def case_id_for(case_dir: Path, repo_root: Path) -> str:
    try:
        rel = case_dir.resolve().relative_to(repo_root.resolve())
    except ValueError:
        rel = Path(case_dir.name)

    parts = rel.parts
    if parts and parts[0] == "example":
        parts = parts[1:]
    return "/".join(parts) or case_dir.name


def case_path_for(case_dir: Path, repo_root: Path) -> str:
    try:
        rel = case_dir.resolve().relative_to(repo_root.resolve())
        path = rel.as_posix()
    except ValueError:
        path = case_dir.as_posix()
    return path.rstrip("/") + "/"


def git_describe(case_dir: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "describe", "--tags", "--always", "--dirty"],
            cwd=case_dir,
            text=True,
            capture_output=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None

    value = result.stdout.strip()
    return value or None


def parse_wallclock_seconds(logfile: Path) -> float | None:
    if not logfile.exists():
        return None

    patterns = (
        re.compile(r"total\s+elapsed\s+time(?:\s+in\s+nekton)?\s*=?\s*([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)", re.I),
        re.compile(r"elapsed\s+time\s*=?\s*([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)", re.I),
        re.compile(r"usertime\s*=?\s*([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)", re.I),
        re.compile(r"\b([0-9]+(?:\.[0-9]*)?)user\b", re.I),
    )

    wallclock = None
    for line in logfile.read_text(errors="replace").splitlines():
        for pattern in patterns:
            match = pattern.search(line)
            if match:
                wallclock = float(match.group(1))
                break
    return wallclock


def last_numeric_value(line: str) -> float | None:
    matches = re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", line)
    if not matches:
        return None
    return float(matches[-1])


def parse_final_residual(case_dir: Path) -> float | None:
    files = sorted({case_dir / "residual.dat", *case_dir.glob("residu*.dat")})
    for path in reversed(files):
        if not path.exists() or not path.is_file():
            continue
        lines = [line.strip() for line in path.read_text(errors="replace").splitlines() if line.strip()]
        if not lines:
            continue
        value = last_numeric_value(lines[-1])
        if value is not None:
            return value
    return None


def parse_parameters(case_dir: Path) -> dict[str, object]:
    par_files = sorted(case_dir.glob("*.par"))
    if not par_files:
        return {}

    parser = configparser.ConfigParser()
    parser.read(par_files[0])
    if not parser.has_section("GENERAL"):
        return {}

    params: dict[str, object] = {}
    general = parser["GENERAL"]
    if "reynoldsnumber" in general:
        try:
            params["reynolds"] = float(general["reynoldsnumber"])
        except ValueError:
            pass

    for key, value in general.items():
        if key.lower().startswith("userparam"):
            params[key.lower()] = value

    return params


def artifact_role(path: Path) -> str:
    name = path.name
    if name.startswith("BF_") or re.search(r"\.f0\d{4}$", name):
        return "CHECKPOINT"
    if re.fullmatch(r"residu.*\.dat", name):
        return "HISTORY"
    if name == "logfile":
        return "LOG"
    if name.endswith(".png"):
        return "FIGURE"
    if name.endswith(".mp4"):
        return "ANIMATION"
    return "OTHER"


def scan_artifacts(case_dir: Path) -> list[dict[str, object]]:
    paths: set[Path] = set()
    for pattern in ARTIFACT_GLOBS:
        paths.update(Path(match) for match in glob.glob(str(case_dir / pattern)))

    artifacts = []
    for path in sorted(paths, key=lambda item: item.name):
        if not path.is_file():
            continue
        artifacts.append(
            {
                "role": artifact_role(path),
                "path": path.name,
                "size_bytes": path.stat().st_size,
                "mtime_utc": utc_mtime(path),
            }
        )
    return artifacts


def build_manifest(case_dir: Path, jobid: int | None, ranks: int | None) -> dict[str, object]:
    repo_root = repo_root_from_script()
    logfile = case_dir / "logfile"

    return {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id_for(case_dir, repo_root),
        "case_path": case_path_for(case_dir, repo_root),
        "run_date_utc": utc_now(),
        "sbatch_jobid": jobid,
        "ranks": ranks,
        "wallclock_seconds": parse_wallclock_seconds(logfile),
        "final_residual": parse_final_residual(case_dir),
        "mks_build_sha": git_describe(case_dir),
        "parameters": parse_parameters(case_dir),
        "artifacts": scan_artifacts(case_dir),
        "slurm_log": "logfile" if logfile.exists() else None,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write a per-case nekStab run provenance manifest.")
    parser.add_argument("case_dir", help="Path to the case directory")
    parser.add_argument("--jobid", type=int, default=None, help="Slurm job id")
    parser.add_argument("--ranks", type=int, default=None, help="MPI rank count")
    parser.add_argument("--out", default=None, help="Output manifest path")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    case_dir = Path(args.case_dir)
    if not case_dir.exists() or not case_dir.is_dir():
        print(f"error: case_dir does not exist: {case_dir}", file=os.sys.stderr)
        return 1

    out_path = Path(args.out) if args.out else case_dir / "manifest.json"
    manifest = build_manifest(case_dir, args.jobid, args.ranks)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
