#!/bin/bash

## GCC

# sudo apt install -y libmpich-dev mpich libopenblas-dev

export FC="mpif77"
export CC="mpicc"
export FFLAGS="-mcmodel=large -march=native -ffixed-line-length-none -w -Wall -g -fbacktrace"
#export FFLAGS="-Og -ffixed-line-length-none -fbacktrace -Wall -Wextra -Warray-bounds -Warray-temporaries -Wconversion -Wno-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ggdb -fsanitize=leak -fsanitize=address"
export USR_LFLAGS+="-L/usr/lib -lopenblas"

source ${NEKSTAB_SOURCE_ROOT}/core/nekStab.sh
