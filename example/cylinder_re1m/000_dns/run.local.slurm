#!/bin/bash
#SBATCH --job-name=cylinder_re1m-dns
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '3\n2.0 0.0 10.0\n5.0 0.0 10.0\n10.0 0.0 10.0\n' > 1cyl.his
mpiexec -np "$SLURM_NTASKS" ./nek5000
