#!/bin/bash
#SBATCH --job-name=cyl-dmd
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
# hpts reads the probe list from the session .his file and writes the probe values after it, in the same file.
# The file is gitignored, so each run writes it here.
printf '1\n5 0 0\n' > ${CASE}.his

# Phase 1 of 2: DNS from the restart file (t = 150) to t = 200. It writes the
# 100 snapshots ${CASE}0.f00001 ... ${CASE}0.f00100 (one every 0.5 time units)
# that the modal analysis of phase 2 reads. ${CASE}.par holds the modal-analysis
# parameters; the DNS parameters are in dns.par.
# The restart file and dns.par are kept in ../600_modal_pod.
for f in rst_1cyl0.f00001 rst_1cyl0.f00001.nek5000 dns.par; do
  [[ -f "${f}" ]] || cp "../600_modal_pod/${f}" .
done
cp "${CASE}.par" modal.par
trap 'cp -f modal.par "${CASE}.par"; rm -f modal.par' EXIT
cp dns.par "${CASE}.par"
mpiexec -np "${NPROCS}" "${EXE}"

# Phase 2 of 2: modal analysis of the snapshots.
cp modal.par "${CASE}.par"
mpiexec -np "${NPROCS}" "${EXE}"

echo "finished=$(date --iso-8601=seconds)"
