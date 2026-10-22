#!/bin/bash
#SBATCH --job-name=movcyl-newton
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=12:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "1cyl" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file. The file is gitignored, so each run writes it here.
printf '1\n5 0 0\n' > 1cyl.his

# The seed is the last field of the forced DNS (30 periods). The field files are not tracked.
if [[ ! -f BF_seed_1cyl0.f00001 ]]; then
  cp ../010_dns_forced/1cyl0.f00030 BF_seed_1cyl0.f00001
fi

mpiexec -np "${SLURM_NTASKS:-8}" ./nek5000
