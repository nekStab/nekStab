#!/usr/bin/env python3
"""submit_case.py — PRE-check and submit a Nek5000 case to local Slurm.

Usage:
    python3 scripts/submit_case.py <case_dir> [--dry-run] [--force]
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
STATE_FILE = REPO_ROOT / "validation" / "campaign_state.json"
INSPECT_SCRIPT = REPO_ROOT / "scripts" / "inspect_nek_field.py"

REQUIRED_FILES = ["SIZE", "NEKSTAB.inc", "makefile_usr.inc", "run.local.slurm", "SESSION.NAME"]


def err(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def parse_par(par_path: Path) -> dict:
    """Parse Nek5000 .par file (INI-like, # comments, case-insensitive keys)."""
    result: dict = {}
    for line in par_path.read_text().splitlines():
        line = line.split("#")[0].strip()
        if not line or line.startswith("["):
            continue
        if "=" in line:
            key, _, val = line.partition("=")
            result[key.strip().lower()] = val.strip()
    return result


def parse_lx1_from_size(size_path: Path) -> int | None:
    """Extract lx1 from SIZE file."""
    text = size_path.read_text()
    m = re.search(r"parameter\s*\(\s*lx1\s*=\s*(\d+)", text, re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


def inspect_field_lx1(field_path: Path) -> int | None:
    """Run inspect_nek_field.py --json and return lx1."""
    result = subprocess.run(
        [sys.executable, str(INSPECT_SCRIPT), "--json", str(field_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    try:
        data = json.loads(result.stdout)
        return int(data["lx1"])
    except (json.JSONDecodeError, KeyError):
        return None


def load_state() -> dict:
    if not STATE_FILE.exists():
        return {
            "schema_version": 1,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "active_jobs": [],
            "completed_jobs": [],
            "failed_jobs": [],
        }
    with STATE_FILE.open("r") as f:
        return json.load(f)


def save_state(state: dict, fh) -> None:
    fh.seek(0)
    fh.truncate()
    json.dump(state, fh, indent=2)
    fh.write("\n")
    fh.flush()


def append_active_job(entry: dict) -> None:
    STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with STATE_FILE.open("a+") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            fh.seek(0)
            content = fh.read()
            state = json.loads(content) if content.strip() else load_state()
            state["active_jobs"].append(entry)
            save_state(state, fh)
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="PRE-check and submit a Nek5000 case to local Slurm.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("case_dir", help="Path to Nek5000 case directory")
    parser.add_argument("--dry-run", action="store_true", help="Print checks + sbatch command; do not submit")
    parser.add_argument("--force", action="store_true", help="Skip executable + IC existence checks")
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    if not case_dir.is_absolute():
        case_dir = Path.cwd() / case_dir
    case_dir = case_dir.resolve()

    # --- PRE CHECK 1: directory exists ---
    if not case_dir.is_dir():
        err(f"Directory not found: {case_dir}")
        sys.exit(1)

    # --- PRE CHECK 2: required files ---
    missing = [f for f in REQUIRED_FILES if not (case_dir / f).exists()]
    if missing:
        err(f"Missing required files in {case_dir}: {', '.join(missing)}")
        sys.exit(1)

    par_files = list(case_dir.glob("*.par"))
    usr_files = list(case_dir.glob("*.usr"))
    if not par_files:
        err(f"No .par file found in {case_dir}")
        sys.exit(1)
    if not usr_files:
        err(f"No .usr file found in {case_dir}")
        sys.exit(1)

    par_path = par_files[0]
    case_name = par_path.stem

    # --- PRE CHECK 3: parse .par ---
    par = parse_par(par_path)
    start_from = par.get("startfrom", "")
    viscosity_str = par.get("viscosity")
    uparam01_str = par.get("userparam01")
    endtime_str = par.get("endtime")

    try:
        viscosity = float(viscosity_str) if viscosity_str is not None else None
    except ValueError:
        viscosity = None
    try:
        uparam01 = float(uparam01_str) if uparam01_str is not None else None
    except ValueError:
        uparam01 = None
    try:
        endtime = float(endtime_str) if endtime_str is not None else None
    except ValueError:
        endtime = None

    par_summary = {
        "startFrom": start_from,
        "viscosity": viscosity,
        "userParam01": uparam01,
        "endTime": endtime,
    }

    # --- PRE CHECK 4: IC file ---
    ic_file: Path | None = None
    if start_from and start_from != "0":
        ic_file = case_dir / start_from
        if not ic_file.exists():
            err(f"IC file not found: {ic_file}\n  (startFrom = {start_from})")
            sys.exit(1)
    else:
        warn(f"startFrom = '{start_from}' → cold start assumed, skipping IC check.")

    # --- PRE CHECK 5: lx1 consistency ---
    if ic_file is not None and not args.force:
        lx1_size = parse_lx1_from_size(case_dir / "SIZE")
        lx1_ic = inspect_field_lx1(ic_file)
        if lx1_size is None:
            warn("Could not parse lx1 from SIZE file; skipping polynomial-order check.")
        elif lx1_ic is None:
            warn(f"Could not read lx1 from {ic_file.name}; skipping polynomial-order check.")
        elif lx1_size != lx1_ic:
            err(
                f"Polynomial-order mismatch: SIZE says lx1={lx1_size} but "
                f"IC file {ic_file.name} has lx1={lx1_ic}. "
                f"Rebuild nek5000 with the correct SIZE or use a matching IC."
            )
            sys.exit(1)

    # --- PRE CHECK 6: executable ---
    nek_exe = case_dir / "nek5000"
    if not args.force:
        if not nek_exe.exists() or not os.access(nek_exe, os.X_OK):
            err(
                f"./nek5000 not found or not executable in {case_dir}.\n"
                f"  Run `mks {case_name}` to build first."
            )
            sys.exit(1)

    # --- PRE CHECK 7: run.local.slurm readable ---
    slurm_script = case_dir / "run.local.slurm"
    if not os.access(slurm_script, os.R_OK):
        err(f"run.local.slurm is not readable: {slurm_script}")
        sys.exit(1)

    # --- SESSION.NAME stamping ---
    session_name_path = case_dir / "SESSION.NAME"
    abs_path_str = str(case_dir) + "/"
    session_content = f"{case_name}\n{abs_path_str}\n"
    session_name_path.write_text(session_content)

    # --- Relative path for state ---
    try:
        rel_case_dir = str(case_dir.relative_to(REPO_ROOT))
    except ValueError:
        rel_case_dir = str(case_dir)

    # --- Summary ---
    print(f"PRE-check summary for: {rel_case_dir}")
    print(f"  case_name   : {case_name}")
    print(f"  startFrom   : {start_from or '(cold start)'}")
    print(f"  viscosity   : {viscosity}")
    print(f"  userParam01 : {uparam01}")
    print(f"  endTime     : {endtime}")
    if ic_file:
        print(f"  IC file     : {ic_file.name} (found)")
    print(f"  SESSION.NAME: stamped ({case_name} / {abs_path_str})")
    print(f"  nek5000     : {'found' if nek_exe.exists() else 'NOT FOUND (--force)'}")

    sbatch_cmd = f"cd {case_dir} && sbatch run.local.slurm"
    print(f"\nsbatch command: {sbatch_cmd}")

    if args.dry_run:
        print("\n[DRY RUN] sbatch not invoked. Exit 0.")
        sys.exit(0)

    # --- SUBMIT ---
    result = subprocess.run(
        ["sbatch", "run.local.slurm"],
        cwd=str(case_dir),
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        err(f"sbatch failed (exit {result.returncode}):\n{result.stderr}")
        sys.exit(1)

    m = re.search(r"Submitted batch job (\d+)", result.stdout)
    if not m:
        err(f"Could not parse job id from sbatch output:\n{result.stdout}")
        sys.exit(1)

    job_id = m.group(1)
    submitted_at = datetime.now(timezone.utc).isoformat()

    # --- Persist ---
    entry = {
        "case_dir": rel_case_dir,
        "case_name": case_name,
        "job_id": job_id,
        "submitted_at": submitted_at,
        "status": "QUEUED",
        "par_summary": par_summary,
    }
    append_active_job(entry)

    print(f"\nSubmitted job {job_id}")
    print(f"  Monitor : squeue -j {job_id}   (or: qs)")
    print(f"  Check   : python3 scripts/check_case.py {rel_case_dir}")


if __name__ == "__main__":
    main()
