"""Tests for evidence badges, stale flags, blockers, and provenance rendering.

"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATION_DIR = REPO_ROOT / "validation"
INDEX_HTML = VALIDATION_DIR / "index.html"


def _build() -> str:
    result = subprocess.run(
        [sys.executable, str(VALIDATION_DIR / "build_gallery.py")],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, (
        f"build_gallery.py failed: {result.stderr}"
    )
    return INDEX_HTML.read_text(encoding="utf-8")


def test_evidence_badges_present() -> None:
    html = _build()
    count = html.count('badge-status')
    assert count > 0, f"expected badge-status badges, found 0"


def test_badge_blocked_present() -> None:
    html = _build()
    count = html.count('badge-blocked')
    assert count > 0, f"expected at least one badge-blocked pill"


def test_badge_stale_class_in_css() -> None:
    html = _build()
    assert '.badge-stale' in html, "CSS class .badge-stale not found in output"


def test_evidence_row_div_present() -> None:
    html = _build()
    count = html.count('class="evidence-row"')
    assert count > 0, f"expected evidence-row divs, found 0"


def test_provenance_row_present() -> None:
    html = _build()
    count = html.count('class="provenance-row"')
    assert count > 0, f"expected provenance-row divs, found 0"


def test_filter_chips_still_present() -> None:
    html = _build()
    assert 'section-nav' in html or 'family-nav' in html, \
        "section-nav / family-nav missing — filter chips broken"


def test_family_nav_still_present() -> None:
    html = _build()
    assert 'id="family-nav"' in html, "family-nav missing — family nav broken"
