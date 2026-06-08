import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "generate_thumbnail.py"

def run(args):
    return subprocess.run([sys.executable, str(SCRIPT)] + args, capture_output=True, text=True)

def test_help():
    r = run(["--help"])
    assert r.returncode == 0
    assert "thumbnail" in r.stdout.lower()

def test_missing_input():
    r = run(["/nonexistent/path.mp4"])
    assert r.returncode == 2

def test_ffmpeg_absent_handled():
    # When ffmpeg not on PATH, script must exit 1 (not crash)
    if shutil.which("ffmpeg") is not None:
        import pytest
        pytest.skip("ffmpeg installed — cannot test absent-path branch")
    # Create a real path so we get past the file-existence check
    tmp = Path("/tmp/dummy_for_test.mp4")
    tmp.write_bytes(b"not a real mp4")
    try:
        r = run([str(tmp)])
        assert r.returncode == 1
        assert "ffmpeg" in (r.stderr + r.stdout).lower()
    finally:
        tmp.unlink(missing_ok=True)
