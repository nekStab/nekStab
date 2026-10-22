#!/bin/bash
#SBATCH --job-name=cylinder_re100-315_stability_direct_meanflow
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=02:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME

# The base flow is the mean flow of the stage 140. The field file is not tracked.
if [[ ! -f BF_1cyl0.f00001 ]]; then
  cp ../140_baseflow_meanflow/BFA1cyl0.f00001 BF_1cyl0.f00001
fi

mpiexec -np "$SLURM_NTASKS" ./nek5000
