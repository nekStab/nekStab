#!/usr/bin/env bash
#SBATCH --job-name=flip_flop_re62-dns
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=2000
#SBATCH --time=06:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "2cyl" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
