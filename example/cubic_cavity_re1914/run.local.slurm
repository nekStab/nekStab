#!/usr/bin/env bash
#SBATCH --job-name=cav1914-nwt
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
# Re=1914 < Re_c (~1916.6) so the target state is a STEADY fixed point:
# userParam01=2.0 Newton-GMRES (k_dim=100) refines the Re=2500 continuation
# seed (BF_cav0.f00001) down to the Re=1914 baseflow. endTime=1 per matvec
# keeps each Jacobian apply short, so this is far cheaper than the Floquet
# stages; 08:00:00 gives convergence headroom. Output writes nwtcav/BF fields
# that the downstream stability stages will start from.
#SBATCH --time=08:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
mpiexec -np "$SLURM_NTASKS" ./nek5000
