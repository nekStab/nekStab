#!/bin/bash
#SBATCH --job-name=re180-animate
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
# Floquet mode animation (userParam01=4.52): advances the UPO base flow one
# period via nek_advance, writing the leading Floquet mode reconstructed at each
# output frame. Single-period post-process, so a short walltime suffices.
#SBATCH --time=04:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
