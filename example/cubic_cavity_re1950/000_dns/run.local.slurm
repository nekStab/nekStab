#!/bin/bash
#SBATCH --job-name=cav-dns
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
#SBATCH --time=06:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n0.25 0.25 0.0\n' > cav.his

# The start field is the on-orbit checkpoint of ../000_dns_period at t = 1700 (field 17 of its run). The field files
# are not tracked: run ../000_dns_period first.
if [[ ! -f ic_sat.f00001 ]]; then
  cp ../000_dns_period/cav0.f00017 ic_sat.f00001
fi

mpiexec -np "$SLURM_NTASKS" ./nek5000
