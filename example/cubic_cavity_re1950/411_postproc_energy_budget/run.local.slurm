#!/bin/bash
#SBATCH --job-name=cav-budget
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4000
# The budget stores one period of the base flow (7516 steps, 2.5 GB per rank on 16 ranks), builds it once, and then
# marches the real and the imaginary part of the mode over it.
#SBATCH --time=04:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME

# The base flow is the on-orbit snapshot of 000_dns_period at t = 1700 with the period in its time stamp, as in the
# Floquet stage. The field files are not tracked: run ../000_dns_period first.
if [[ ! -f BF_cav0.f00001 ]]; then
  python3 ../000_dns_period/stamp_period.py ../000_dns_period/cav0.f00017 10.78909 BF_cav0.f00001
fi

# The mode and its exponent come from the Floquet stage: run ../311_stability_direct_floquet first.
for f in dRecav0.f00001 dImcav0.f00001 Spectre_NSd_conv.dat; do
  if [[ ! -f "$f" ]]; then
    cp "../311_stability_direct_floquet/$f" .
  fi
done

mpiexec -np "$SLURM_NTASKS" ./nek5000
