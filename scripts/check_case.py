#!/usr/bin/env python3
"""check_case.py -- Poll Slurm and run post-processing for a submitted Nek5000 case.

Usage:
    python3 scripts/check_case.py <case_dir> [--plots] [--strict]
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
STATE_FILE = REPO_ROOT / "validation" / "campaign_state.json"
INSPECT_SCRIPT = REPO_ROOT / "scripts" / "inspect_nek_field.py"
PLOT_RUNNER = REPO_ROOT / "scripts" / "plot_runner.py"

PLOT_SCRIPTS = [
    "p_residu.py",
    "p_residu_nwt.py",
    "p_spec.py",
    "plot.py",
    "p_lift.py",
    "p_his.py",
    "p_energy.py",
]

ERROR_MARKERS = ["error", "fatal", "diverged", "nan"]
COMPLETION_MARKERS = ["total elapsed", "step      "]


def err(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)


def load_state_locked(fh) -> dict:
    fh.seek(0)
    content = fh.read()
    if content.strip():
        return json.loads(content)
    return {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "active_jobs": [],
        "completed_jobs": [],
        "failed_jobs": [],
    }


def save_state_locked(state: dict, fh) -> None:
    fh.seek(0)
    fh.truncate()
    json.dump(state, fh, indent=2)
    fh.write("\n")
    fh.flush()


def fmt_elapsed(delta_seconds: float) -> str:
    h = int(delta_seconds // 3600)
    m = int((delta_seconds % 3600) // 60)
    s = int(delta_seconds % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def check_slurm(job_id: str, submitted_at: str) -> str | None:
    """Return None if job is finished, else a status string."""
    try:
        result = subprocess.run(
            ["squeue", "-j", job_id, "-h", "-o", "%T"],
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return None

    stdout = result.stdout.strip()
    stderr = result.stderr.strip()

    if not stdout or "Invalid job id" in stderr or "Invalid job id" in stdout:
        return None

    if stdout == "PD":
        return "PENDING"
    if stdout == "R":
        try:
            sub_dt = datetime.fromisoformat(submitted_at)
            elapsed = (datetime.now(timezone.utc) - sub_dt).total_seconds()
            return f"RUNNING  (elapsed: {fmt_elapsed(elapsed)})"
        except Exception:
            return "RUNNING"
    return f"{stdout}  (still in flight)"


def detect_outcome(logfile: Path) -> tuple[str, str, str]:
    """Return (outcome, convergence_reason, log_tail)."""
    if not logfile.exists() or logfile.stat().st_size == 0:
        return "ambiguous", "logfile missing or empty", ""

    try:
        lines = logfile.read_text(errors="replace").splitlines()
    except Exception as exc:
        return "ambiguous", f"logfile unreadable: {exc}", ""

    tail_lines = lines[-50:]
    log_tail = "\n".join(tail_lines)[:2000]
    lower_tail = "\n".join(tail_lines).lower()

    for marker in ERROR_MARKERS:
        if marker in lower_tail:
            for line in tail_lines:
                if marker in line.lower():
                    return "failed", line.strip(), log_tail

    for marker in COMPLETION_MARKERS:
        if marker in lower_tail:
            return "completed", "No error markers found", log_tail

    return "completed", "No error markers found", log_tail


def get_new_field_files(case_dir: Path, after_ts: float) -> list[Path]:
    new_files = []
    for f in sorted(case_dir.glob("*.f0????")):
        if f.stat().st_mtime > after_ts:
            new_files.append(f)
    return new_files


def run_inspect(field_path: Path) -> str:
    result = subprocess.run(
        [sys.executable, str(INSPECT_SCRIPT), str(field_path)],
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def run_plot_script(script_path: Path, case_dir: Path) -> None:
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    if script_path.name == "plot.py":
        subprocess.run(
            [sys.executable, str(PLOT_RUNNER), str(case_dir), "--timeout", "300"],
            cwd=str(case_dir),
            env=env,
            timeout=330,
        )
    else:
        print(f"  Running {script_path.name} ...")
        try:
            result = subprocess.run(
                [sys.executable, script_path.name],
                cwd=str(case_dir),
                env=env,
                capture_output=True,
                text=True,
                timeout=300,
            )
            if result.returncode != 0:
                print(f"    FAILED (exit {result.returncode}): {result.stderr[:200]}")
            else:
                print("    OK")
        except subprocess.TimeoutExpired:
            print("    TIMEOUT after 300s")


def normalize_case_dir(case_dir_str: str) -> str:
    p = Path(case_dir_str)
    if not p.is_absolute():
        p = Path.cwd() / p
    p = p.resolve()
    try:
        return str(p.relative_to(REPO_ROOT))
    except ValueError:
        return str(p)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Poll Slurm and run post-processing for a submitted Nek5000 case.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("case_dir", help="Path to Nek5000 case directory")
    parser.add_argument("--plots", action="store_true", help="Run plot scripts found in case dir")
    parser.add_argument("--strict", action="store_true", help="Exit non-zero on FAILED or AMBIGUOUS")
    args = parser.parse_args()

    rel_case_dir = normalize_case_dir(args.case_dir)
    abs_case_dir = (REPO_ROOT / rel_case_dir).resolve()

    if not STATE_FILE.exists():
        err(f"State file not found: {STATE_FILE}")
        sys.exit(1)

    with STATE_FILE.open("r") as f:
        state = json.load(f)

    matches = [j for j in state.get("active_jobs", []) if j.get("case_dir") == rel_case_dir]
    if not matches:
        err(f"No active job for {rel_case_dir}")
        sys.exit(1)

    entry = sorted(matches, key=lambda j: j.get("submitted_at", ""))[-1]
    job_id = entry["job_id"]
    submitted_at = entry["submitted_at"]

    print(f"Checking case: {rel_case_dir}  (job {job_id})")

    status = check_slurm(job_id, submitted_at)
    if status is not None:
        print(f"  Slurm status: {status}")
        sys.exit(0)

    print(f"  Job {job_id} is no longer in the Slurm queue -- checking outcome.")

    logfile = abs_case_dir / "logfile"
    outcome, convergence_reason, log_tail = detect_outcome(logfile)
    print(f"  Outcome: {outcome.upper()}")
    print(f"  Reason : {convergence_reason}")

    try:
        sub_dt = datetime.fromisoformat(submitted_at)
        after_ts = sub_dt.timestamp()
    except Exception:
        after_ts = 0.0

    new_field_files = get_new_field_files(abs_case_dir, after_ts)
    print(f"  New field files ({len(new_field_files)}):")
    for ff in new_field_files:
        print(f"    {ff.name}")
        info = run_inspect(ff)
        for line in info.splitlines():
            print(f"      {line}")

    new_plot_files: list[Path] = []
    if args.plots:
        print("  Running plot scripts:")
        png_before = {f: f.stat().st_mtime for f in abs_case_dir.glob("*.png")}
        for script_name in PLOT_SCRIPTS:
            sp = abs_case_dir / script_name
            if sp.exists():
                run_plot_script(sp, abs_case_dir)
        for f in sorted(abs_case_dir.glob("*.png")):
            old_mtime = png_before.get(f, 0.0)
            if f.stat().st_mtime > old_mtime:
                new_plot_files.append(f)
        print(f"  New/updated plots: {[f.name for f in new_plot_files]}")

    with STATE_FILE.open("r+") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            state2 = load_state_locked(fh)
            matches2 = [j for j in state2.get("active_jobs", []) if j.get("case_dir") == rel_case_dir]
            if matches2:
                to_move = sorted(matches2, key=lambda j: j.get("submitted_at", ""))[-1]
                state2["active_jobs"] = [j for j in state2["active_jobs"] if j is not to_move]
                if to_move in state2["active_jobs"]:
                    state2["active_jobs"].remove(to_move)
            else:
                to_move = dict(entry)

            to_move["finished_at"] = datetime.now(timezone.utc).isoformat()
            to_move["outcome"] = outcome
            to_move["new_field_files"] = [str(f.relative_to(abs_case_dir)) for f in new_field_files]
            to_move["new_plot_files"] = [str(f.relative_to(abs_case_dir)) for f in new_plot_files]
            to_move["convergence_reason"] = convergence_reason
            to_move["last_log_tail"] = log_tail

            dest_key = "failed_jobs" if outcome in ("failed", "ambiguous") else "completed_jobs"
            state2.setdefault(dest_key, []).append(to_move)
            save_state_locked(state2, fh)
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)

    print(f"  State updated: moved to {dest_key}.")

    if args.strict and outcome in ("failed", "ambiguous"):
        sys.exit(1)


if __name__ == "__main__":
    main()
