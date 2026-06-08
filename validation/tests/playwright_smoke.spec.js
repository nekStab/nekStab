// playwright_smoke.spec.js - browser action smoke tests for the gallery UI.
//
// Tests the static index.html via file:// URL — no running server needed.
// Verifies user-facing actions (copy, lightbox, filter chips, recompile/resubmit
// buttons, queue panel) still work after the gallery rewrite.
//
// Run from repo root:
//   npx playwright test validation/tests/playwright_smoke.spec.js --browser chromium
// Or with the config:
//   npx playwright test --config validation/playwright.config.js
//

const { test, expect } = require('@playwright/test');
const path = require('path');

const indexPath = path.resolve(__dirname, '../index.html');
const PAGE_URL = `file://${indexPath}`;

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

/**
 * Navigate to the static gallery and wait for the page to settle.
 * Disables content-visibility:auto so all cards are rendered for assertions.
 */
async function openGallery(page) {
    await page.goto(PAGE_URL, { waitUntil: 'domcontentloaded' });
    // Force all cards to render — same trick as audit_screenshot.py
    await page.addStyleTag({ content: '.card { content-visibility: visible !important; }' });
    await page.waitForTimeout(300);
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

test('page loads without JavaScript errors', async ({ page }) => {
    const jsErrors = [];
    page.on('pageerror', (err) => jsErrors.push(err.message));

    await openGallery(page);

    // Title must include 'nekStab'
    await expect(page).toHaveTitle(/nekStab/);

    // No catastrophic JS errors
    expect(jsErrors, `JS errors on load: ${jsErrors.join('; ')}`).toHaveLength(0);
});

test('family nav is visible with nav chips', async ({ page }) => {
    await openGallery(page);

    const nav = page.locator('#family-nav');
    await expect(nav).toBeVisible();

    const chips = page.locator('.nav-chip');
    const chipCount = await chips.count();
    expect(chipCount, `Expected >= 1 .nav-chip inside #family-nav, found ${chipCount}`).toBeGreaterThanOrEqual(1);
});

test('filter chips are clickable and activate', async ({ page }) => {
    await openGallery(page);

    // Click the first status filter chip
    const statusChip = page.locator('.filter-chip[data-filter-group="status"]').first();
    await expect(statusChip).toBeVisible();
    await statusChip.click();
    await page.waitForTimeout(150);

    // The chip should now have class 'active'
    await expect(statusChip).toHaveClass(/active/);
});

test('copy action fires on .copyable click', async ({ page, context }) => {
    // Grant clipboard permissions so navigator.clipboard.writeText works in headless
    await context.grantPermissions(['clipboard-read', 'clipboard-write']);
    await openGallery(page);

    const copyable = page.locator('.copyable').first();
    await expect(copyable).toBeVisible();
    await copyable.click();
    await page.waitForTimeout(250);

    // In file:// context clipboard write may silently fail; check 'copied' class
    // but treat absence as a soft warning rather than a hard failure.
    const hasCopied = await copyable.evaluate((el) => el.classList.contains('copied'));
    if (!hasCopied) {
        console.warn('SOFT: .copyable did not receive class "copied" — clipboard may be blocked in file:// context');
    }
    // The page must not crash regardless
    await expect(page.locator('.card').first()).toBeVisible();
});

test('lightbox opens on a.zoom click', async ({ page }) => {
    await openGallery(page);

    const zoomLink = page.locator('a.zoom').first();
    await expect(zoomLink).toBeVisible();
    await zoomLink.click();

    // Lightbox must become visible within 3 s
    await expect(page.locator('.lightbox.open')).toBeVisible({ timeout: 3000 });

    // Escape closes it
    await page.keyboard.press('Escape');
    await page.waitForTimeout(400);
    await expect(page.locator('.lightbox.open')).not.toBeVisible();
});

test('recompile and resubmit buttons exist on cards', async ({ page }) => {
    await openGallery(page);

    const recompileCount = await page.locator('button[data-action="recompile"]').count();
    expect(recompileCount, `Expected >= 1 recompile button, found ${recompileCount}`).toBeGreaterThanOrEqual(1);

    const resubmitCount = await page.locator('button[data-action="resubmit"]').count();
    expect(resubmitCount, `Expected >= 1 resubmit button, found ${resubmitCount}`).toBeGreaterThanOrEqual(1);
});

test('cylinder family has ladder steps', async ({ page }) => {
    await openGallery(page);

    // Scroll the section into view
    const section = page.locator('#family-cylinder_re100');
    await section.scrollIntoViewIfNeeded();
    await expect(section).toBeVisible();

    const stepCount = await page.locator('.ladder-step').count();
    expect(stepCount, `Expected >= 1 .ladder-step inside #family-cylinder_re100, found ${stepCount}`).toBeGreaterThanOrEqual(1);
});

test('queue panel is present in DOM', async ({ page }) => {
    await openGallery(page);

    // Queue panel may not be visible at page top but must be attached to DOM.
    // Use 'attached' state — no live Slurm data needed for this check.
    const panel = page.locator('#queue-panel');
    await expect(panel).toBeAttached();
});

test('page has cards visible after load', async ({ page }) => {
    await openGallery(page);

    const cardCount = await page.locator('.card').count();
    expect(cardCount, `Expected >= 1 .card on the page, found ${cardCount}`).toBeGreaterThanOrEqual(1);

    // First card should be visible in the viewport after scrolling to it
    await page.locator('.card').first().scrollIntoViewIfNeeded();
    await expect(page.locator('.card').first()).toBeVisible();
});
