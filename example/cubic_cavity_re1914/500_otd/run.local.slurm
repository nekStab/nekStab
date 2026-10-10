#!/bin/bash
#SBATCH --job-name=cubic_cavity_re1914-500_otd
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
#SBATCH --time=08:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n0.25 0.25 0.0\n' > cav.his

# The OTD basis evolves about the steady base flow of 210_baseflow_newton/fp, which the DNS does
# not advance. The field file is not tracked.
if [[ ! -f BF_cav0.f00001 ]]; then
  cp ../210_baseflow_newton/fp/BF_cav0.f00001 .
fi

# The OTD driver reads r01cav0.f<n> and r02cav0.f<n> (the highest n found) as the initial
# perturbations, and a run writes files with these names. Remove those of an earlier run, or
# the driver starts from them. The seeds are the real and imaginary parts of the leading
# eigenmode of 310_stability_direct, so the basis starts on the invariant subspace of the
# leading pair. OTD_SEED=0 starts from white noise instead.
rm -f r[0-9][0-9]cav0.f*
if [[ "${OTD_SEED:-1}" == "1" ]]; then
  cp ../310_stability_direct/dRecav0.f00001 r01cav0.f00001
  cp ../310_stability_direct/dImcav0.f00001 r02cav0.f00001
fi

mpiexec -np "$SLURM_NTASKS" ./nek5000
