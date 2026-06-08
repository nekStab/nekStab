#!/usr/bin/env python3
"""audit_screenshot.py - drive the gallery via Playwright for UX audit evidence.

Captures:
  shots/desktop.png       - full page, 1440x900
  shots/desktop-hd.png    - full page, 1920x1080
  shots/mobile.png        - full page, 390x844 (iPhone-ish)
  shots/hover.png         - one card hovered
  shots/lightbox.png      - lightbox open on first thumbnail
  shots/copy.png          - one .copyable in 'copied' state
  shots/console.log       - browser console messages
  shots/network.log       - request/response log
  shots/family-nav.png    - family navigation bar
  shots/filter-bar.png    - filter chip bar
  shots/cylinder-family.png - Static Cylinder family section (first 800px)

Run via: uv run --with playwright python validation/audit_screenshot.py
Requires `python -m playwright install chromium` once (done).
Server must already be on http://127.0.0.1:8000/ (validation/serve.py).
"""
from __future__ import annotations
from pathlib import Path
from playwright.sync_api import sync_playwright, Page  # pyright: ignore[reportMissingImports]

OUT = Path(__file__).resolve().parent / 'shots'
OUT.mkdir(exist_ok=True)
URL = 'http://127.0.0.1:8000/'


def assert_present(page: Page, selector: str, label: str):
    """Assert that at least one visible element matching *selector* exists.

    Raises AssertionError with a clear 'MISSING: ...' message so CI logs
    immediately identify the broken selector rather than a cryptic timeout.
    """
    locator = page.locator(selector)
    count = locator.count()
    if count == 0:
        message = f'MISSING: {label} ({selector})'
        print(message)
        raise AssertionError(message)
    if not locator.first.is_visible():
        message = f'MISSING: {label} ({selector}) present but not visible'
        print(message)
        raise AssertionError(message)
    print(f'  found {count}x {selector} ✓')
    return locator


def main() -> None:
    console_log: list[str] = []
    network_log: list[str] = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)

        ctx = browser.new_context(viewport={'width': 1440, 'height': 900})
        page = ctx.new_page()

        page.on('console', lambda m: console_log.append(
            f'[{m.type}] {m.text}  (from {m.location.get("url","?")}:{m.location.get("lineNumber","?")})'
        ))
        page.on('requestfailed', lambda r: network_log.append(
            f'FAIL {r.method} {r.url} -- {r.failure}'
        ))
        page.on('response', lambda r: network_log.append(
            f'{r.status} {r.request.method} {r.url}'
        ))

        # 1. desktop 1440
        page.goto(URL, wait_until='networkidle')
        # disable content-visibility:auto for full-page captures (it skips
        # rendering offscreen cards, which screenshots without scroll).
        page.add_style_tag(content='.card { content-visibility: visible !important; }')
        page.wait_for_timeout(400)  # let entrance animations settle

        # --- selector assertions: fail early with clear messages ---
        print('checking selectors...')
        assert_present(page, '#family-nav', 'family nav')
        assert_present(page, '#family-nav .nav-chip', 'nav chips in #family-nav')
        assert_present(page, '#filter-bar', 'filter bar')
        assert_present(page, '#filter-bar .filter-chip', 'filter chips')
        assert_present(page, '.badge-status', 'evidence badges')
        assert_present(page, '#family-cylinder_re100', 'Static Cylinder family')
        assert_present(page, '#family-cylinder_re100 .ladder-step', 'cylinder_re100 ladder steps')
        print('all selectors OK')

        page.screenshot(path=str(OUT / 'desktop.png'), full_page=True)
        print('desktop ✓', OUT / 'desktop.png')

        # 2. desktop hd 1920
        page.set_viewport_size({'width': 1920, 'height': 1080})
        page.wait_for_timeout(200)
        page.screenshot(path=str(OUT / 'desktop-hd.png'), full_page=True)
        print('desktop-hd ✓', OUT / 'desktop-hd.png')

        # 3. mobile
        page.set_viewport_size({'width': 390, 'height': 844})
        page.wait_for_timeout(200)
        page.screenshot(path=str(OUT / 'mobile.png'), full_page=True)
        print('mobile ✓', OUT / 'mobile.png')

        # 4. hover on a card
        page.set_viewport_size({'width': 1440, 'height': 900})
        page.wait_for_timeout(200)
        first_card = page.locator('.card').first
        first_card.hover()
        page.wait_for_timeout(400)
        page.screenshot(path=str(OUT / 'hover.png'), clip={'x': 0, 'y': 280, 'width': 1440, 'height': 620})
        print('hover ✓', OUT / 'hover.png')

        # 5. lightbox open
        # click an actual <a.zoom> (not just .card — some cards have no PNG)
        zoom_link = page.locator('a.zoom').first
        zoom_link.click()
        page.wait_for_selector('.lightbox.open', timeout=2000)
        page.wait_for_timeout(400)  # fade-in settle
        page.screenshot(path=str(OUT / 'lightbox.png'), full_page=False)
        print('lightbox ✓', OUT / 'lightbox.png')
        page.keyboard.press('Escape')
        page.wait_for_timeout(300)

        # 6. copy state on a .copyable
        copyable = page.locator('.copyable').first
        # grant clipboard so navigator.clipboard.writeText works in headless
        ctx.grant_permissions(['clipboard-read', 'clipboard-write'])
        copyable.click()
        page.wait_for_timeout(100)
        page.screenshot(path=str(OUT / 'copy.png'), clip={'x': 0, 'y': 280, 'width': 1440, 'height': 360})
        print('copy ✓', OUT / 'copy.png')

        # 7. queue panel close-up
        queue = page.locator('#queue-panel')
        queue.scroll_into_view_if_needed()
        bbox = queue.bounding_box()
        if bbox:
            page.screenshot(path=str(OUT / 'queue-panel.png'), clip={
                'x': max(0, bbox['x'] - 8),
                'y': max(0, bbox['y'] - 8),
                'width': min(1440, bbox['width'] + 16),
                'height': bbox['height'] + 16,
            })
            print('queue-panel ✓', OUT / 'queue-panel.png')

        # 8. family nav bar close-up
        # scroll back to top and capture the sticky nav
        page.evaluate('window.scrollTo(0, 0)')
        page.wait_for_timeout(200)
        nav = page.locator('#family-nav')
        nav.scroll_into_view_if_needed()
        nav_bbox = nav.bounding_box()
        if nav_bbox:
            page.screenshot(path=str(OUT / 'family-nav.png'), clip={
                'x': max(0, nav_bbox['x'] - 8),
                'y': max(0, nav_bbox['y'] - 8),
                'width': min(1440, nav_bbox['width'] + 16),
                'height': nav_bbox['height'] + 16,
            })
            print('family-nav ✓', OUT / 'family-nav.png')

        # 9. filter bar close-up
        fbar = page.locator('#filter-bar')
        fbar.scroll_into_view_if_needed()
        fbar_bbox = fbar.bounding_box()
        if fbar_bbox:
            page.screenshot(path=str(OUT / 'filter-bar.png'), clip={
                'x': max(0, fbar_bbox['x'] - 8),
                'y': max(0, fbar_bbox['y'] - 8),
                'width': min(1440, fbar_bbox['width'] + 16),
                'height': fbar_bbox['height'] + 16,
            })
            print('filter-bar ✓', OUT / 'filter-bar.png')

        # 10. cylinder family section (first 800 px of height — captures ladder)
        cyl = page.locator('#family-cylinder_re100')
        cyl.scroll_into_view_if_needed()
        page.wait_for_timeout(300)
        cyl_bbox = cyl.bounding_box()
        if cyl_bbox:
            page.screenshot(path=str(OUT / 'cylinder-family.png'), clip={
                'x': max(0, cyl_bbox['x'] - 8),
                'y': max(0, cyl_bbox['y'] - 8),
                'width': min(1440, cyl_bbox['width'] + 16),
                'height': min(800, cyl_bbox['height'] + 16),
            })
            print('cylinder-family ✓', OUT / 'cylinder-family.png')

        browser.close()

    (OUT / 'console.log').write_text('\n'.join(console_log) + '\n')
    (OUT / 'network.log').write_text('\n'.join(network_log) + '\n')
    print(f'\nconsole: {len(console_log)} msgs')
    print(f'network: {len(network_log)} events')
    fails = [n for n in network_log if n.startswith('FAIL') or n.startswith('4') or n.startswith('5')]
    if fails:
        print(f'\n⚠ {len(fails)} failed/error responses:')
        for f in fails[:10]:
            print('  ', f)


if __name__ == '__main__':
    main()
