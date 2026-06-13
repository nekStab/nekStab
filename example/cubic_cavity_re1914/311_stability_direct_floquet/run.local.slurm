#!/usr/bin/env bash
#SBATCH --job-name=cubic_cavity_re1914-311_stability_direct_floquet
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
