#!/bin/bash
#SBATCH --job-name=cav-floq
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
# One Floquet matvec integrates the full period (~7496 steps, ~270 s on 16
# ranks); a k_dim=192 Arnoldi factorization alone is ~14.6 h, so the previous
# 06:00:00 limit killed the run mid-factorization at vector 92/192 (#68). Allow
# a full factorization plus restart headroom.
#SBATCH --time=24:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
