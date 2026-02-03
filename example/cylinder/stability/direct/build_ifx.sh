#!/bin/bash
# Build with ifx (no fp-model=consistent)
set -e

cd ~/nekStab/example/cylinder/stability/direct

# Clean
rm -rf obj nek5000 drive.o makefile makefile_usr.inc NEKSTAB.inc .state .usr 2>/dev/null || true
rm -rf ~/nekStab/Nek5000/3rd_party/gslib/lib 2>/dev/null || true

# Source Intel
source /opt/intel/oneapi/setvars.sh --force

# Set ifx
export NEKSTAB_FC=ifx
echo "Using: $(which ifx)"
echo "Version: $(ifx --version | head -1)"

# Build
mks 1cyl
