#!/bin/bash
#SBATCH --job-name=re180-dfloq
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
# Direct Floquet (userParam01=3.11): the monodromy map is built by integrating
# the converged UPO base flow (BF_1cyl0.f00001) over one shedding period per
# matvec; k_dim=100 (uparam07) Arnoldi vectors -> allow generous walltime.
#SBATCH --time=12:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
