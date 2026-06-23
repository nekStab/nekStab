#!/bin/bash
#SBATCH --job-name=naca-310dir
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=04:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "naca0012" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
