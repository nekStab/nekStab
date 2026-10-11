#!/bin/bash
set -uo pipefail
# Rebuild the two RANS stability variants, for example after a change of src/eigensolvers.f90. Run it from any directory.
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../../.." && pwd)"
cd "$root"
[ -f bin/sourceme.sh ] && source bin/sourceme.sh 2>/dev/null || true
export PATH="$root/bin:$PATH"
for d in coupled quasilaminar; do
  cd "$here/$d"
  rm -f obj/eigensolvers.o obj/libnek5000.a nek5000
  echo "=== building $d ==="
  if mks 1cyl > build_eig.log 2>&1 && [ -x nek5000 ]; then
    echo "$d: BUILD OK  ($(ls -la nek5000 | awk '{print $5, $6, $7, $8}'))"
    grep -c "eigensolvers" build_eig.log | sed 's/^/  eigensolvers recompiled lines: /'
  else
    echo "$d: BUILD FAIL"; grep -iE "error|undefined" build_eig.log | tail -15
  fi
done
echo "=== DONE ==="
