#!/bin/bash

## GCC

# sudo apt install -y libmpich-dev mpich libopenblas-dev

export FC="mpif77"
export CC="mpicc"
export FFLAGS="-mcmodel=large -march=native -ffixed-line-length-none -g -fbacktrace"
export USR_LFLAGS+="-L/usr/lib -lopenblas"

source ${NEKSTAB_SOURCE_ROOT}/core/nekStab.sh
