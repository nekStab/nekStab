#!/usr/bin/env bash
# check_session_name.sh — verify or normalize SESSION.NAME in each case dir.
#
# Each Nek5000 case dir must contain a SESSION.NAME of the form:
#     <basename>
#     <absolute path to case dir>/
#
# where <basename> matches the local .usr / .par / .re2 / .ma2 files.
#
# Usage:
#   scripts/check_session_name.sh                  # audit only (default)
#   scripts/check_session_name.sh --fix            # rewrite mismatched files
#
# Exit codes:
#   0  every case OK
#   1  one or more mismatches (audit) or rewrites (fix)
#   2  usage error
set -euo pipefail

MODE="audit"
case "${1:-}" in
    ""|--audit|-a)     MODE="audit" ;;
    --fix|-f)          MODE="fix" ;;
    -h|--help)
        sed -n '2,18p' "$0"; exit 0 ;;
    *)
        echo "usage: $0 [--audit|--fix]" >&2; exit 2 ;;
esac

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FAIL=0
TOUCH=0

while IFS= read -r -d '' sn; do
    dir="$(dirname "$sn")"
    abs_dir="$(cd "$dir" && pwd)/"

    # Detect basename: prefer .re2 mesh (unique per case); fall back to .par.
    re2="$(find "$dir" -maxdepth 1 -name '*.re2' -not -path '*/obj/*' | head -1)"
    if [[ -n "$re2" ]]; then
        base="$(basename "$re2" .re2)"
    else
        par="$(find "$dir" -maxdepth 1 -name '*.par' -not -path '*/obj/*' | head -1)"
        if [[ -z "$par" ]]; then
            echo "SKIP  $sn  (no .re2 or .par sibling)"
            continue
        fi
        base="$(basename "$par" .par)"
    fi

    expected=$'%s\n%s\n'
    # shellcheck disable=SC2059
    expected="$(printf "$expected" "$base" "$abs_dir")"
    actual="$(cat "$sn")"

    if [[ "$actual" == "$expected" ]]; then
        echo "OK    $sn"
    else
        FAIL=$((FAIL+1))
        echo "DIFF  $sn"
        if [[ "$MODE" == "fix" ]]; then
            printf '%s\n%s\n' "$base" "$abs_dir" > "$sn"
            TOUCH=$((TOUCH+1))
            echo "FIX   $sn  ->  $base | $abs_dir"
        fi
    fi
done < <(find "$ROOT/example" -name 'SESSION.NAME' -not -path '*/_templates/*' -print0)

echo
if [[ "$MODE" == "fix" ]]; then
    echo "Summary: $FAIL diffs, $TOUCH rewritten."
    exit $(( TOUCH > 0 ? 1 : 0 ))
else
    echo "Summary: $FAIL mismatch(es)."
    exit $(( FAIL > 0 ? 1 : 0 ))
fi
