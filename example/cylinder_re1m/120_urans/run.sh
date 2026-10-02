#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
mpiexec -np 8 ./nek5000 > logfile 2>&1
