"""The gallery draws every stage of example/ with its figures, and nothing else.

validation/build_gallery.py builds validation/index.html from the stage folders. These tests
fail when

- a stage folder has no reference, summary, description, figure or README (one test per stage);
- a stage folder is missing on the page, or the page shows a card without a stage folder;
- the page has a deferred, blocked or pending card;
- a link or an image on the page does not resolve, or points to an untracked file;
- a folder of example/ has no title in build_gallery.FAMILIES;
- the mode names of the gallery differ from the mode codes of src/mode_codes.f90;
- validation/index.html differs from the output of build_gallery.py.
"""

from __future__ import annotations

import re
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
VALIDATION = REPO / "validation"
INDEX = VALIDATION / "index.html"
sys.path.insert(0, str(VALIDATION))

import build_gallery as gallery  # noqa: E402

STAGES = gallery.discover_stages()


def tracked_files() -> set[str] | None:
    """Return the files that git tracks, or None outside a git checkout."""
    if not (REPO / ".git").exists():
        return None
    listing = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "-z"], capture_output=True, check=True
    ).stdout.decode()
    return set(listing.split("\0"))


class PageParser(HTMLParser):
    """Collect the references, the ids and the cards of the page. Text of style and script is dropped."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.references: list[tuple[str, str]] = []
        self.ids: set[str] = set()
        self.cards: list[dict[str, str]] = []
        self.thumbs: list[dict[str, str]] = []
        self.text: list[str] = []
        self._skip = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = {k: v or "" for k, v in attrs}
        if tag in ("style", "script"):
            self._skip += 1
        if "id" in values:
            self.ids.add(values["id"])
        for key in ("src", "href", "data-img-src"):
            if values.get(key):
                self.references.append((tag, values[key]))
        classes = values.get("class", "").split()
        if tag == "div" and "card" in classes:
            self.cards.append(values)
        if tag == "button" and "thumb" in classes:
            self.thumbs.append(values)

    def handle_endtag(self, tag: str) -> None:
        if tag in ("style", "script"):
            self._skip -= 1

    def handle_data(self, data: str) -> None:
        if not self._skip:
            self.text.append(data)


@pytest.fixture(scope="module")
def page() -> PageParser:
    parser = PageParser()
    parser.feed(INDEX.read_text(encoding="utf-8"))
    return parser


@pytest.mark.parametrize("stage", STAGES, ids=lambda s: s.path)
def test_stage_has_reference_summary_figure_and_readme(stage: gallery.Stage) -> None:
    assert not stage.problems, "; ".join(stage.problems)


def test_there_are_stages() -> None:
    assert len(STAGES) >= 50


def test_every_stage_has_a_card_with_all_its_figures(page: PageParser) -> None:
    cards = {c["data-case-id"] for c in page.cards}
    thumbs = {(t["data-case-id"], t["data-img-src"]) for t in page.thumbs}
    for stage in gallery.drawn_stages(STAGES):
        assert stage.path in cards, f"{stage.path} has no card"
        for figure in stage.figures:
            assert (stage.path, gallery._img_url(figure)) in thumbs, (
                f"{figure} is not on the page"
            )


def test_every_card_has_a_stage_folder(page: PageParser) -> None:
    folders = {s.path for s in STAGES}
    assert not [
        c["data-case-id"] for c in page.cards if c["data-case-id"] not in folders
    ]


def test_no_deferred_blocked_or_pending_card(page: PageParser) -> None:
    html = INDEX.read_text(encoding="utf-8")
    assert not re.search(r'class="card[^"]*\b(deferred|blocked|pending)\b', html)
    assert not re.search(r"data-status=", html)
    text = " ".join(page.text).lower()
    for word in ("deferred", "blocked", "pending", "not yet", "todo"):
        assert word not in text, f"'{word}' on the page"
    assert page.cards, "no card on the page"


def test_no_dead_link(page: PageParser) -> None:
    tracked = tracked_files()
    for tag, reference in page.references:
        if "://" in reference or reference.startswith(("mailto:", "/")):
            continue
        if reference.startswith("#"):
            assert reference[1:] in page.ids, (
                f"<{tag}> points to the missing id {reference}"
            )
            continue
        target = (VALIDATION / reference.split("#", 1)[0].split("?", 1)[0]).resolve()
        assert target.is_file(), f"<{tag}> points to the missing file {reference}"
        if tracked is not None:
            assert target.relative_to(REPO).as_posix() in tracked, (
                f"{reference} is not tracked by git"
            )


def test_every_folder_of_example_has_a_title() -> None:
    folders = {s.family for s in STAGES}
    assert not sorted(folders - set(gallery.FAMILIES)), (
        "folders without a title in FAMILIES"
    )
    assert not sorted(set(gallery.FAMILIES) - folders), (
        "titles in FAMILIES without a stage folder"
    )


def test_mode_names_match_the_mode_codes_of_the_solver() -> None:
    text = (REPO / "src" / "mode_codes.f90").read_text()
    codes = {int(c) for c in re.findall(r"::\s*MODE_\w+\s*=\s*(\d+)", text)}
    names = {round(float(key) * 100) for key in gallery.MODE_NAMES}
    assert names == codes, (
        f"only in the gallery: {sorted(names - codes)}; only in the solver: {sorted(codes - names)}"
    )


def test_every_mode_of_a_stage_has_a_name() -> None:
    unknown = {
        s.path: s.mode for s in STAGES if s.mode and s.mode not in gallery.MODE_NAMES
    }
    assert not unknown


def test_the_committed_page_is_the_output_of_the_generator() -> None:
    assert INDEX.read_text(encoding="utf-8") == gallery.build_page(), (
        "validation/index.html is out of date: run python validation/build_gallery.py"
    )
