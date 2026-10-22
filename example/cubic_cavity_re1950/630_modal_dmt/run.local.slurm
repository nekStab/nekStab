#!/bin/bash
#SBATCH --job-name=cav-dmt
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
#SBATCH --time=04:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME

# The start field is the on-orbit checkpoint of 000_dns_period at t = 1700. The field files are not tracked.
if [[ ! -f ic_sat.f00001 ]]; then
  if [[ -f ../000_dns/ic_sat.f00001 ]]; then
    cp ../000_dns/ic_sat.f00001 .
  else
    cp ../000_dns_period/cav0.f00017 ic_sat.f00001
  fi
fi

mpiexec -np "$SLURM_NTASKS" ./nek5000
