#!/bin/bash
#SBATCH --job-name=cubic_cavity_re1920-210_baseflow_newton
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3000
# Re=1920 > Re_c (1916.63): the steady flow is unstable, so a DNS cannot reach it and Newton-GMRES
# (userParam01=2.0, k_dim=100) must find it. The seed BF_cav0.f00001 is the last field of ../../000_dns (Re=1900).
# endTime=1 per Jacobian application keeps each matvec short. Newton writes the converged flow over the
# same file name; it is the base flow of ../../310_stability_direct.
#SBATCH --time=08:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "cav" "$(pwd)" > SESSION.NAME
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n0.25 0.25 0.0\n' > cav.his
mpiexec -np "$SLURM_NTASKS" ./nek5000
