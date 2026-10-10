"""Tests for the parts of a card of the gallery: reference badge and table, provenance line, navigation."""

from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATION_DIR = REPO_ROOT / "validation"
sys.path.insert(0, str(VALIDATION_DIR))

import build_gallery as gallery  # noqa: E402


@lru_cache(maxsize=1)
def _build() -> str:
    return gallery.build_page()


def test_reference_badge_on_every_card() -> None:
    html = _build()
    assert html.count('class="badge-status') == html.count('<div class="card"')
    assert html.count('class="badge-status') > 0


def test_reference_table_on_every_card() -> None:
    html = _build()
    assert html.count('<details class="ref">') == html.count('<div class="card"')
    assert 'class="ref-table"' in html


def test_evidence_row_div_present() -> None:
    assert 'class="evidence-row"' in _build()


def test_provenance_row_present() -> None:
    assert 'class="provenance-row"' in _build()


def test_filter_chips_still_present() -> None:
    html = _build()
    assert 'id="filter-bar"' in html
    assert 'data-filter-group="lane"' in html


def test_family_nav_still_present() -> None:
    assert 'id="family-nav"' in _build()


def test_build_is_deterministic() -> None:
    assert gallery.build_page() == gallery.build_page()


def test_quantity_values_are_formatted() -> None:
    assert gallery._format_number(0.1247875) == "0.1247875"
    assert gallery._format_number(8.496308e-06) == "8.496308e-06"
    assert gallery._format_number(12) == "12"
    assert gallery._tolerance_text({"tol_rel": 0.05}) == "±0.05 rel"
    assert gallery._tolerance_text({"tol_abs": 1e-3}) == "±0.001"
