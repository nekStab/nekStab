"""plot_runner.py — Wrap plot.py invocations with logging, timeout, and env snapshot.

Usage:
    python3 scripts/plot_runner.py <case_dir> [--timeout 600] [--log plot.log] [--env-snapshot]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def _write_env_snapshot(case_dir: Path) -> None:
    """Write Python/matplotlib/numpy versions and pip freeze to env.log."""
    env_log = case_dir / "env.log"
    lines: list[str] = []

    # Python version
    lines.append(f"python: {sys.version}")

    # matplotlib version
    try:
        import importlib.metadata
        lines.append(f"matplotlib: {importlib.metadata.version('matplotlib')}")
    except Exception:
        lines.append("matplotlib: (not importable)")

    # numpy version
    try:
        import importlib.metadata
        lines.append(f"numpy: {importlib.metadata.version('numpy')}")
    except Exception:
        lines.append("numpy: (not importable)")

    # pip freeze
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        lines.append("")
        lines.append("--- pip freeze ---")
        lines.append(result.stdout)
    except Exception as exc:
        lines.append(f"pip freeze failed: {exc}")

    env_log.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run plot.py inside a case directory with logging and timeout.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("case_dir", help="Path to directory containing plot.py")
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        metavar="SECONDS",
        help="Subprocess timeout in seconds (default: 600)",
    )
    parser.add_argument(
        "--log",
        default=None,
        metavar="FILENAME",
        help="Log file path (default: <case_dir>/plot.log)",
    )
    parser.add_argument(
        "--env-snapshot",
        action="store_true",
        help="Write Python/matplotlib/numpy versions and pip freeze to <case_dir>/env.log",
    )
    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    if not case_dir.is_dir():
        print(f"ERROR: case_dir does not exist: {case_dir}", file=sys.stderr)
        sys.exit(1)

    log_path = Path(args.log) if args.log else case_dir / "plot.log"

    if args.env_snapshot:
        _write_env_snapshot(case_dir)

    # Build subprocess environment with MPLBACKEND=Agg
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"

    t0 = time.monotonic()
    exit_code = 0
    status_word = "PASS"

    try:
        with log_path.open("w") as log_fh:
            result = subprocess.run(
                [sys.executable, "plot.py"],
                cwd=case_dir,
                env=env,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
            )
        elapsed = time.monotonic() - t0
        if result.returncode != 0:
            status_word = "FAIL"
            exit_code = 1
        else:
            status_word = "PASS"
            exit_code = 0
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        status_word = "TIMEOUT"
        exit_code = 2
        with log_path.open("a") as log_fh:
            log_fh.write(f"\n[plot_runner] TIMEOUT after {args.timeout}s\n")
    except Exception as exc:
        elapsed = time.monotonic() - t0
        status_word = "FAIL"
        exit_code = 1
        with log_path.open("a") as log_fh:
            log_fh.write(f"\n[plot_runner] ERROR: {exc}\n")

    print(f"{status_word} [{args.case_dir}] in {elapsed:.1f}s")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
