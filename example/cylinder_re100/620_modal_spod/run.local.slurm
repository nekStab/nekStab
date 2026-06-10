#!/usr/bin/env bash
#SBATCH --job-name=cyl-spod
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2000
#SBATCH --time=04:00:00
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail

CASE=${CASE:-1cyl}
NPROCS=${SLURM_NTASKS:-4}
EXE=${EXE:-./nek5000}

# Slurm opens --output with truncation by default. Outside Slurm, re-exec
# through tee so manual runs of this script also use the same local log name.
if [[ -z "${SLURM_JOB_ID:-}" && -z "${LOGFILE_REDIRECTED:-}" ]]; then
  export LOGFILE_REDIRECTED=1
  exec > >(tee logfile) 2>&1
fi

echo "host=$(hostname)"
echo "case=${CASE}"
echo "ntasks=${NPROCS}"
echo "workdir=$(pwd)"
echo "started=$(date --iso-8601=seconds)"

if [[ "${NEK_DRY_RUN:-0}" == "1" ]]; then
  echo "dry_run=1"
  mpiexec -np "${NPROCS}" hostname
  echo "finished=$(date --iso-8601=seconds)"
  exit 0
fi

if [[ ! -x "${EXE}" ]]; then
  echo "ERROR: executable not found or not executable: ${EXE}" >&2
  echo "Build first with: mks ${CASE}" >&2
  exit 1
fi

if [[ ! -f SESSION.NAME ]]; then
  echo "ERROR: SESSION.NAME not found in $(pwd)" >&2
  exit 1
fi

if [[ ! -f "${CASE}.par" ]]; then
  echo "ERROR: ${CASE}.par not found in $(pwd)" >&2
  exit 1
fi

printf "%s\n%s/\n" "${CASE}" "$(pwd)" > SESSION.NAME

mpiexec -np "${NPROCS}" "${EXE}"

echo "finished=$(date --iso-8601=seconds)"
