#!/bin/bash
# Minimal script to run Nek5000 case and run restart logic
# Emulates supercomputer PBS script for local use

# Load environment
# conda activate nek

# Compile the case (uncomment if you want to always recompile)
# mks 1cyl --fresh

echo 1cyl > SESSION.NAME
echo $PWD'/' >> SESSION.NAME

# Run the simulation (replace 1cyl and 5 as needed)
mpiexec -n 5 ./nek5000 > logfile 2> logerror

# Post-processing: adjust/restart logic
python3 check_restart.py > logrestart 2>&1
