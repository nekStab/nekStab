// playwright.config.js — minimal Playwright config for the gallery smoke tests.
//
// Run from repo root:
//   ./validation/node_modules/.bin/playwright test --config validation/playwright.config.js
//
// Requires: npm install in validation/ (installs @playwright/test locally).
// DO NOT install @playwright/test globally — keep it local to validation/.

module.exports = {
    testDir: './tests',
    testMatch: ['**/playwright_smoke.spec.js'],
    use: {
        headless: true,
        browserName: 'chromium',
    },
    // No retries — smoke tests should pass or fail cleanly.
    retries: 0,
    // One worker avoids flaky parallel clipboard races.
    workers: 1,
};
