#!/bin/bash
#SBATCH --job-name=movcyl-dns
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=24:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
CASE=${CASE:-1cyl}
NPROCS=${SLURM_NTASKS:-8}
EXE=${EXE:-./nek5000}
if [[ -z "${SLURM_JOB_ID:-}" && -z "${LOGFILE_REDIRECTED:-}" ]]; then
  export LOGFILE_REDIRECTED=1
  exec > >(tee logfile) 2>&1
fi
printf "%s\n%s/\n" "${CASE}" "$(pwd)" > SESSION.NAME
mpiexec -np "${NPROCS}" "${EXE}"
