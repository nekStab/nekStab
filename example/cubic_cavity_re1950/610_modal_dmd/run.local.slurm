#!/bin/bash
#SBATCH --job-name=cav-dmd
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

# The 100 snapshots are the fields 80 to 179 of ../000_dns (t = 1740 to 1789.5, one every 0.5 time units).
# The field files are not tracked: run ../000_dns first.
for i in $(seq 1 100); do
  link=$(printf 'cav0.f%05d' "$i")
  [[ -e "$link" ]] || ln -sf "../000_dns/$(printf 'cav0.f%05d' $((i + 79)))" "$link"
done

mpiexec -np "$SLURM_NTASKS" ./nek5000
