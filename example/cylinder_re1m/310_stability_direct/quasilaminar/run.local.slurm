#!/bin/bash
#SBATCH --job-name=re1m-quasilaminar
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n5 0 0\n' > 1cyl.his
mpiexec -np "$SLURM_NTASKS" ./nek5000
