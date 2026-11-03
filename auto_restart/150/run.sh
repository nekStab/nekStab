#!/bin/bash
# Run one job of the restart chain: Nek5000 until the endTime of the .par file,
# then check_restart.py, which backs up the outputs and starts the next job.
# Emulates a supercomputer batch script for local use.
#
# Needs NEKSTAB_SOURCE_ROOT in the environment (check_restart.py imports
# bin/nekStab_tools.py from it) and an executable ./nek5000 built with mks.
# NP is the number of MPI ranks (at least lpmin of SIZE); change it in the
# environment or here.

# Load environment
# conda activate nek

# Compile the case (uncomment if you want to always recompile)
# mks <case> --fresh

# The case name is the name of the .usr file.
CASE=$(basename "$(ls ./*.usr | head -1)" .usr)
NP=${NP:-5}

echo "$CASE" > SESSION.NAME
echo $PWD'/' >> SESSION.NAME

# Run the simulation
mpiexec -n "$NP" ./nek5000 > logfile 2> logerror
# Tell check_restart.py that the solver has ended, and how: a crashed first job
# leaves a field file of time 0, which looks like a first start without this file.
echo $? > logstatus

# Post-processing: adjust/restart logic
python3 check_restart.py > logrestart 2>&1
