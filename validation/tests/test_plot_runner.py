"""Tests for scripts/plot_runner.py."""
from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

PLOT_RUNNER = Path(__file__).resolve().parent.parent.parent / "scripts" / "plot_runner.py"


def _make_plot_py(case_dir: Path, content: str) -> None:
    (case_dir / "plot.py").write_text(textwrap.dedent(content))


def test_plot_runner_success(tmp_path: Path) -> None:
    _make_plot_py(tmp_path, """\
        print("hello from plot")
        import sys; sys.exit(0)
    """)
    result = subprocess.run(
        [sys.executable, str(PLOT_RUNNER), str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    assert "PASS" in result.stdout
    log_text = (tmp_path / "plot.log").read_text()
    assert "hello from plot" in log_text


def test_plot_runner_timeout(tmp_path: Path) -> None:
    _make_plot_py(tmp_path, """\
        import time
        time.sleep(30)
    """)
    result = subprocess.run(
        [sys.executable, str(PLOT_RUNNER), str(tmp_path), "--timeout", "1"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    assert "TIMEOUT" in result.stdout


def test_plot_runner_failure(tmp_path: Path) -> None:
    _make_plot_py(tmp_path, """\
        import sys; sys.exit(1)
    """)
    result = subprocess.run(
        [sys.executable, str(PLOT_RUNNER), str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    assert "FAIL" in result.stdout


def test_env_snapshot(tmp_path: Path) -> None:
    _make_plot_py(tmp_path, "pass\n")
    result = subprocess.run(
        [sys.executable, str(PLOT_RUNNER), str(tmp_path), "--env-snapshot"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"
    env_log = tmp_path / "env.log"
    assert env_log.exists(), "env.log not created"
    env_text = env_log.read_text()
    assert "matplotlib" in env_text, f"matplotlib not in env.log: {env_text!r}"


def test_custom_log_path(tmp_path: Path) -> None:
    _make_plot_py(tmp_path, 'print("custom log test")\n')
    custom_log = tmp_path / "custom.log"
    result = subprocess.run(
        [sys.executable, str(PLOT_RUNNER), str(tmp_path), "--log", str(custom_log)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert custom_log.exists()
    assert "custom log test" in custom_log.read_text()
