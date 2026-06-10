#!/usr/bin/env bash
#SBATCH --job-name=cyl150-dns
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

if [[ -z "${SLURM_JOB_ID:-}" && -z "${LOGFILE_REDIRECTED:-}" ]]; then
  export LOGFILE_REDIRECTED=1
  exec > >(tee logfile) 2>&1
fi

echo "host=$(hostname)"; echo "case=${CASE}"; echo "ntasks=${NPROCS}"
echo "workdir=$(pwd)"; echo "started=$(date --iso-8601=seconds)"

[[ -x "${EXE}" ]] || { echo "ERROR: missing executable" >&2; exit 1; }
[[ -f SESSION.NAME ]] || { echo "ERROR: SESSION.NAME missing" >&2; exit 1; }
[[ -f "${CASE}.par" ]] || { echo "ERROR: ${CASE}.par missing" >&2; exit 1; }

printf "%s\n%s/\n" "${CASE}" "$(pwd)" > SESSION.NAME
mpiexec -np "${NPROCS}" "${EXE}"
echo "finished=$(date --iso-8601=seconds)"
