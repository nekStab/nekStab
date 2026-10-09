#!/bin/bash
#SBATCH --job-name=re180-otd
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=06:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile
set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n5 0 0\n' > 1cyl.his

# The DNS starts on the orbit of 210_baseflow_newton (the field file is not tracked).
if [[ ! -f BF_1cyl0.f00001 ]]; then
  cp ../210_baseflow_newton/BF_1cyl0.f00001 .
fi

# The OTD driver reads r01<case>0.f<n> (the highest n found) as the initial perturbation and
# starts from white noise when there is none. A run writes files with these names, so remove
# those of an earlier run, or the driver starts from them.
rm -f r[0-9][0-9]1cyl0.f*
mpiexec -np "$SLURM_NTASKS" ./nek5000
