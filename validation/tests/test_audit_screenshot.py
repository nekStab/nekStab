"""Static selector tests for validation/audit_screenshot.py.

These tests intentionally avoid importing Playwright or requiring a running
validation server. They compile the audit script and inspect validation/index.html
directly so selector regressions fail in plain pytest.
"""

from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path
import py_compile

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATION_DIR = REPO_ROOT / "validation"
INDEX_HTML = VALIDATION_DIR / "index.html"


try:
    from bs4 import BeautifulSoup  # type: ignore[import-untyped]
except ImportError:  # pragma: no cover - exercised only without bs4 installed
    BeautifulSoup = None


class Element:
    def __init__(
        self,
        tag: str,
        attrs: dict[str, str],
        parent: "Element | None" = None,
    ) -> None:
        self.tag = tag
        self.attrs = attrs
        self.parent = parent
        self.children: list[Element] = []

    def get(self, name: str) -> str | None:
        return self.attrs.get(name)

    def has_class(self, class_name: str) -> bool:
        return class_name in (self.attrs.get("class") or "").split()

    def find(self, selector: str | None = None, **kwargs: str) -> "Element | None":
        if selector is None and "id" in kwargs:
            selector = f"#{kwargs['id']}"
        if selector is None:
            return None
        return next(iter(self.select(selector)), None)

    def select(self, selector: str) -> list["Element"]:
        parts = selector.split()
        current = [self]
        for part in parts:
            matches: list[Element] = []
            for node in current:
                matches.extend(
                    child
                    for child in node.iter_descendants()
                    if _matches_selector(child, part)
                )
            current = matches
        return current

    def iter_descendants(self) -> list["Element"]:
        nodes: list[Element] = []
        stack = list(self.children)
        while stack:
            node = stack.pop(0)
            nodes.append(node)
            stack[0:0] = node.children
        return nodes


class FallbackParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.root = Element("[document]", {})
        self.stack = [self.root]

    def handle_starttag(
        self,
        tag: str,
        attrs: list[tuple[str, str | None]],
    ) -> None:
        node = Element(
            tag,
            {name: value or "" for name, value in attrs},
            parent=self.stack[-1],
        )
        self.stack[-1].children.append(node)
        if tag not in {
            "area",
            "base",
            "br",
            "col",
            "embed",
            "hr",
            "img",
            "input",
            "link",
            "meta",
            "param",
            "source",
            "track",
            "wbr",
        }:
            self.stack.append(node)

    def handle_endtag(self, tag: str) -> None:
        for index in range(len(self.stack) - 1, 0, -1):
            if self.stack[index].tag == tag:
                del self.stack[index:]
                return


def _matches_selector(node: Element, selector: str) -> bool:
    if selector.startswith("#"):
        return node.get("id") == selector[1:]
    if selector.startswith("."):
        return node.has_class(selector[1:])
    if "." in selector:
        tag, class_name = selector.split(".", 1)
        return node.tag == tag and node.has_class(class_name)
    return node.tag == selector


def soup():
    html = INDEX_HTML.read_text(encoding="utf-8")
    if BeautifulSoup is not None:
        return BeautifulSoup(html, "html.parser")
    parser = FallbackParser()
    parser.feed(html)
    return parser.root


def test_audit_script_compiles() -> None:
    py_compile.compile(str(VALIDATION_DIR / "audit_screenshot.py"), doraise=True)


def test_index_has_family_nav() -> None:
    assert soup().find(id="family-nav") is not None


def test_index_has_filter_bar() -> None:
    doc = soup()
    filter_bar = doc.find(id="filter-bar")
    assert filter_bar is not None
    assert filter_bar.select(".filter-chip")


def test_index_has_cylinder_ladder() -> None:
    doc = soup()
    cylinder = doc.find(id="family-cylinder_re100")
    assert cylinder is not None
    assert cylinder.select(".ladder-step")


def test_index_has_badge_status() -> None:
    assert soup().select(".badge-status")


def test_index_has_cards() -> None:
    assert soup().select(".card")


def test_index_has_zoom_links() -> None:
    assert soup().select("a.zoom")


def test_index_has_copyable() -> None:
    assert soup().select(".copyable")


def test_shots_dir_structure() -> None:
    assert (VALIDATION_DIR / "shots").is_dir()


def test_family_nav_has_chips() -> None:
    doc = soup()
    family_nav = doc.find(id="family-nav")
    assert family_nav is not None
    assert family_nav.select(".nav-chip")
