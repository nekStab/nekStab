#!/bin/bash
#SBATCH --job-name=tsyphon_bf-D5
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=2000
#SBATCH --time=01:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "tsyphon" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
