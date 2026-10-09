#!/bin/bash
#SBATCH --job-name=cav-dns-unstable
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
# DNS of 2000 time units from the orbit plus a small part of the leading Floquet
# mode: about 1.4e6 steps at dt = 1.44e-3, measured at 725 time units in 5 h on
# 16 ranks, so about 14 h for the whole run. Allow headroom.
#SBATCH --time=24:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
