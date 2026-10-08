#!/bin/bash
#SBATCH --job-name=cav-floq-adj
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
# One Floquet matvec integrates one period (7500 steps, about 6 min on 16 ranks).
# With k_dim = 100 the run took 8.4 h on 16 ranks; k_dim = 192 would not fit in 24 h.
#SBATCH --time=24:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
