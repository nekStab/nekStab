#!/bin/bash
set -e

## GCC

# sudo apt install -y libmpich-dev mpich libopenblas-dev

export FC="mpif77"
export CC="mpicc"
export FFLAGS="-mcmodel=large -march=native -ffixed-line-length-none -w -Wall -g -fbacktrace"
#export FFLAGS="-mcmodel=large -march=native -ffixed-line-length-none -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan"

export USR_LFLAGS+="-L/usr/lib -lopenblas"

source ${NEKSTAB_SOURCE_ROOT}/core/nekStab.sh
