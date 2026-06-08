#!/usr/bin/env bash
# check_mp4_sizes.sh — reject *.mp4 staged for commit if any exceeds 5 MB.
#
# Designed to be called from .git/hooks/pre-commit.
# Usage: check_mp4_sizes.sh [path/to/file.mp4]
#   Without argument: checks all staged .mp4 files.
#   With argument: checks the specified file (useful for smoke-testing).
set -euo pipefail

MAX_SIZE=5242880  # 5 MB in bytes

failed=0

if [ "$#" -ge 1 ]; then
    # Direct file argument (smoke-test mode)
    files="$1"
else
    # Read staged mp4 files from git index
    files=$(git -C "$(dirname "$0")/../.." diff --cached --name-only --diff-filter=AM 2>/dev/null \
            | grep -i '\.mp4$' || true)
fi

if [ -z "$files" ]; then
    exit 0
fi

while IFS= read -r f; do
    [ -z "$f" ] && continue
    if [ ! -f "$f" ]; then
        echo "check_mp4_sizes: WARNING: '$f' not found, skipping." >&2
        continue
    fi
    size=$(stat -c '%s' "$f" 2>/dev/null || stat -f '%z' "$f")
    if [ "$size" -gt "$MAX_SIZE" ]; then
        echo "check_mp4_sizes: FAIL: $f is ${size} bytes (limit: ${MAX_SIZE} bytes = 5 MB)" >&2
        failed=1
    fi
done <<< "$files"

if [ "$failed" -ne 0 ]; then
    echo "check_mp4_sizes: Reduce file size (increase CRF or lower resolution) then re-stage." >&2
    exit 1
fi

exit 0
