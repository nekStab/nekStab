#!/bin/bash

## INTEL COMPILER
export FC="mpiifort"
export CC="mpiicc"
export FFLAGS="-mcmodel=medium -shared-intel -xHost -extend-source -mkl"

# OPTIONAL FLAG
export USR_LFLAGS+="-I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl" #optional link for old versions

source ${NEKSTAB_SOURCE_ROOT}/core/nekStab.sh
