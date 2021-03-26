#!/bin/bash

export PGI=/opt/pgi;
export PATH=/opt/pgi/linux86-64/19.10/bin:$PATH;
export MANPATH=$MANPATH:/opt/pgi/linux86-64/19.10/man;
export LM_LICENSE_FILE=$LM_LICENSE_FILE:/opt/pgi/license.dat; 

## PGI COMPILER
export FC="pgfortran"
export CC="pgcc"
export FFLAGS="-mcmodel=medium -g -traceback"
export USR_LFLAGS+="-llapack -lblas"

source ${NEKSTAB_SOURCE_ROOT}/core/nekStab.sh
