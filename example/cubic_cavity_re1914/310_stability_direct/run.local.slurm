#!/bin/bash
#SBATCH --job-name=cubic_cavity_re1914-310_stability_direct
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n0.25 0.25 0.0\n' > cav.his
mpiexec -np "$SLURM_NTASKS" ./nek5000
