"""Baseline regression test for validation/build_gallery.py.

Locks down the current gallery behavior so the upcoming catalog rewrite
has a known reference.

Asserts:

* validation/build_gallery.py, validation/serve.py, and
  validation/audit_screenshot.py compile (`py_compile`);
* `python3 validation/build_gallery.py` exits 0 and writes
  validation/index.html;
* the rendered HTML contains at least one validated card, at least one
  deferred card, at least one figure image, and at least one
  recompile/resubmit action button.

The assertions check behavior (class/attribute presence), not whitespace,
so deliberate layout tweaks should not break this test.

Runs as plain `python3 validation/tests/test_baseline_gallery.py` or via
pytest.
"""

from __future__ import annotations

import py_compile
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATION_DIR = REPO_ROOT / "validation"
INDEX_HTML = VALIDATION_DIR / "index.html"

SCRIPTS_TO_COMPILE = [
    VALIDATION_DIR / "build_gallery.py",
    VALIDATION_DIR / "serve.py",
    VALIDATION_DIR / "audit_screenshot.py",
]


def test_scripts_compile() -> None:
    for path in SCRIPTS_TO_COMPILE:
        assert path.is_file(), f"missing: {path}"
        py_compile.compile(str(path), doraise=True)


def test_gallery_build_runs() -> None:
    result = subprocess.run(
        [sys.executable, str(VALIDATION_DIR / "build_gallery.py")],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, (
        f"build_gallery.py exited {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert INDEX_HTML.is_file(), "validation/index.html not produced"


def test_gallery_baseline_assertions() -> None:
    html = INDEX_HTML.read_text(encoding="utf-8")

    validated_cards = re.findall(r'<div class="card"(?:\s|>)', html)
    deferred_cards = re.findall(r'<div class="card deferred"(?:\s|>)', html)
    figure_images = re.findall(r'<img [^>]*src="figures/[^"]+"', html)
    action_buttons = re.findall(
        r'data-action="(?:recompile|resubmit)"', html
    )

    assert validated_cards, "no validated cards found"
    assert deferred_cards, "no deferred cards found"
    assert figure_images, "no figure images found"
    assert action_buttons, "no recompile/resubmit action buttons found"


def main() -> int:
    failures: list[str] = []
    for fn in (
        test_scripts_compile,
        test_gallery_build_runs,
        test_gallery_baseline_assertions,
    ):
        try:
            fn()
            print(f"  ok  {fn.__name__}")
        except AssertionError as exc:
            print(f"  FAIL {fn.__name__}: {exc}")
            failures.append(fn.__name__)
        except Exception as exc:
            print(f"  ERROR {fn.__name__}: {exc}")
            failures.append(fn.__name__)
    if failures:
        print(f"\n{len(failures)} test(s) failed: {failures}")
        return 1
    print("\nall baseline gallery checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
