#!/bin/bash
#SBATCH --job-name=movcyl-floquet
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=06:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file. The file is gitignored, so each run writes it here.
printf '1\n5 0 0\n' > 1cyl.his

# The periodic orbit is the output of the Newton stage. The field files are not tracked.
if [[ ! -f BF_1cyl0.f00001 ]]; then
  cp ../210_baseflow_newton/BF_1cyl0.f00001 .
fi

mpiexec -np "${SLURM_NTASKS:-8}" ./nek5000
