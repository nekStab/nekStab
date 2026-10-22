#!/bin/bash
#SBATCH --job-name=cylinder_re100-140_baseflow_meanflow
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=02:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n5 0 0\n' > 1cyl.his

# The start field is the t = 100 field of the DNS stage. The field file is tracked there.
if [[ ! -f rst_1cyl0.f00001 ]]; then
  cp ../000_dns/rst_1cyl0.f00001 .
fi

mpiexec -np "$SLURM_NTASKS" ./nek5000
