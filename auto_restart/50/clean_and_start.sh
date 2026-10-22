#!/bin/bash

ps aux | egrep 'nek5000|python|mpiexec|mpirun' | grep -v grep

# Double-strategy: pgrep+kill, then pkill -f as backup
for proc in nek5000 python python3 mpiexec mpirun; do
    # First, kill all matching PIDs for this user
    pids=$(pgrep -u "$USER" $proc)
    if [ -n "$pids" ]; then
        echo "[INFO] Killing $proc PIDs: $pids"
        kill -9 $pids
    else
        echo "[INFO] No $proc processes found for user $USER via pgrep"
    fi
    # Second, fallback to pkill -9 -f for any stragglers
    pkill -9 -f $proc && echo "[INFO] pkill -9 -f $proc: DONE (backup)" || echo "[INFO] pkill -9 -f $proc: NONE FOUND (backup)"
done

# Wait a moment for processes to terminate
sleep 2
ps aux | egrep 'nek5000|python|mpiexec|mpirun' | grep -v grep

# Clean previous outputs
rm 1cyl0.f00* 1cyl.his_* lift_drag.* log*

# Update .par file for restart
python3 -c "from check_restart import prepare_par_for_restart, parameter_glob, single_job_time, case_preserving_adjust_par_file; prepare_par_for_restart(parameter_glob, single_job_time, case_preserving_adjust_par_file)"

mks 1cyl

# Start the simulation
./run_nohup.sh