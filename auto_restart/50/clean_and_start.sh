#!/bin/bash
# Start a restart chain from scratch in this folder: stop what runs here, delete
# the outputs of the previous chain, point the .par file at rst_<case>0.f00001
# and the first end time, rebuild the case and start job 1.
#
# Needs NEKSTAB_SOURCE_ROOT, NEK_SOURCE_ROOT and mks in the environment.

# The case name is the name of the .usr file.
CASE=$(basename "$(ls ./*.usr | head -1)" .usr)
HERE=$(pwd -P)

# Stop only the processes of this chain: those whose working directory is this
# folder. Killing by name would also kill the user's other solvers and Python
# jobs. The shell running this script and its parent are kept.
chain_pids() {
    # The patterns are anchored so that an editor open on run.sh is not matched.
    for pid in $( { pgrep -u "$USER" -x 'nek5000|mpiexec|mpiexec\.hydra|mpirun|hydra_pmi_proxy'
                    pgrep -u "$USER" -f '^(/usr/bin/|/bin/)?bash (\./)?(run|run_nohup)\.sh'
                    pgrep -u "$USER" -f '^(/usr/bin/)?python3? (\./)?check_restart\.py'
                  } ); do
        [ "$pid" = "$$" ] && continue
        [ "$pid" = "$PPID" ] && continue
        [ "$(readlink "/proc/$pid/cwd" 2>/dev/null)" = "$HERE" ] && echo "$pid"
    done
}

# Cancel the Slurm jobs of this folder, queued or running. A queued job has no
# process yet, so chain_pids does not see it; it would start after the cleanup below
# and write into the same files as the job started at the end of this script.
if command -v squeue > /dev/null; then
    while IFS='|' read -r jobid jobdir; do
        [ -n "$jobid" ] && [ "$(realpath "$jobdir" 2>/dev/null)" = "$HERE" ] && { echo "[INFO] Cancelling Slurm job $jobid of this folder"; scancel "$jobid"; }
    done < <(squeue -h -u "$USER" -o '%i|%Z')
fi

pids=$(chain_pids)
if [ -n "$pids" ]; then
    echo "[INFO] Killing PIDs of this chain: $pids"
    kill -9 $pids
else
    echo "[INFO] No process of this chain found"
fi

# Wait a moment for processes to terminate
sleep 2

# Clean previous outputs (-f: not all of them exist at the first start)
rm -f "${CASE}"0.f00* "${CASE}".his_* lift_drag.* total_energy.* total_enstrophy.* log*

# Update .par file for restart
python3 -c "from check_restart import prepare_par_for_restart, parameter_glob, single_job_time, case_preserving_adjust_par_file; prepare_par_for_restart(parameter_glob, single_job_time, case_preserving_adjust_par_file)"

mks "$CASE"

# Start the simulation: on a Slurm machine through job.pbs, as every later job of the
# chain is started by check_restart.py; otherwise in the background of this shell.
if command -v sbatch > /dev/null; then
    sbatch job.pbs
else
    ./run_nohup.sh
fi
