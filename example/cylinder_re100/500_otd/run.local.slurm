#!/bin/bash
#SBATCH --job-name=cyl-otd
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

if [[ ! -f "${CASE}.par" ]]; then
  echo "ERROR: ${CASE}.par not found in $(pwd)" >&2
  exit 1
fi

printf "%s\n%s/\n" "${CASE}" "$(pwd)" > SESSION.NAME

# The OTD basis evolves about the steady base flow of 310_stability_direct. The
# field file is not in this folder in a fresh clone.
if [[ ! -f "BF_${CASE}0.f00001" ]]; then
  cp "../310_stability_direct/direct/BF_${CASE}0.f00001" .
fi

# The OTD driver reads r01<case>0.f<n> and r02<case>0.f<n> (the highest n found) as the
# initial perturbations, and a run writes files with these names. Remove those of an
# earlier run, or the driver starts from them. The seeds are the real and imaginary
# parts of the leading eigenmode of 310_stability_direct, so the basis starts on the
# leading invariant subspace. OTD_SEED=0 starts from white noise instead (set
# endTime = 200 in the .par file).
rm -f r[0-9][0-9]"${CASE}"0.f*
if [[ "${OTD_SEED:-1}" == "1" ]]; then
  cp "../310_stability_direct/direct/dRe${CASE}0.f00001" "r01${CASE}0.f00001"
  cp "../310_stability_direct/direct/dIm${CASE}0.f00001" "r02${CASE}0.f00001"
fi

mpiexec -np "${NPROCS}" "${EXE}"

echo "finished=$(date --iso-8601=seconds)"
