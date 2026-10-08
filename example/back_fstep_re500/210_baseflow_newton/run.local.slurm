#!/bin/bash
#SBATCH --job-name=bfs-newton
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=02:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "bfs" "$(pwd)" > SESSION.NAME

# Phase 1 of 2: DNS from the uniform initial condition settles the steady
# recirculation on this mesh (writes bfs0.f00001..4).
cp bfs_dns.par bfs.par
mpiexec -np "$SLURM_NTASKS" ./nek5000

# Phase 2 of 2: Newton-GMRES polishes the settled field (startFrom = bfs0.f00004)
# to the steady base flow BF_bfs0.f00001.
cp bfs_newton.par bfs.par
mpiexec -np "$SLURM_NTASKS" ./nek5000

# Leave the directory in its starting state: bfs.par is the phase-1 file.
cp bfs_dns.par bfs.par
