#!/usr/bin/env bash
#SBATCH --job-name=cav1914-nwtfp
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
# Re=1914 < Re_c (~1916.6): the target is a STEADY fixed point. Newton-GMRES
# (userParam01=2.0, k_dim=100) refines the Re=2500 continuation seed
# (BF_cav0.f00001) down to the Re=1914 baseflow. endTime=1 per Jacobian apply
# keeps each matvec short -> far cheaper than the Floquet stages. The converged
# BF_cav0.f00001 written here is the initial condition for the 310/320/330/411/
# modal stages of this geometry.
#SBATCH --time=08:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
