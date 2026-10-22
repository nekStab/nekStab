#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '3\n2.0 0.0 10.0\n5.0 0.0 10.0\n10.0 0.0 10.0\n' > 1cyl.his
mpiexec -np 8 ./nek5000 > logfile 2>&1
