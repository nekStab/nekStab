#!/bin/bash
set -uo pipefail
cd /home/rfrantz/nekStab
[ -f bin/sourceme.sh ] && source bin/sourceme.sh 2>/dev/null || true
export PATH="/home/rfrantz/nekStab/bin:$PATH"
for d in 310_stability_direct/coupled 310_stability_direct/quasilaminar; do
  cd "/home/rfrantz/nekStab/example/cylinder_re1m/$d"
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
