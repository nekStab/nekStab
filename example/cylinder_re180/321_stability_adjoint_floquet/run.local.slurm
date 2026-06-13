#!/usr/bin/env bash
#SBATCH --job-name=re180-afloq
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
# Adjoint Floquet (userParam01=3.21): adjoint monodromy map over one period of
# the UPO base flow (BF_1cyl0.f00001), replayed time-reversed; k_dim=100
# (uparam07) -> allow generous walltime.
#SBATCH --time=12:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
