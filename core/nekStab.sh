#!/bin/bash

#add to .bashrc
#export NEKSTAB_SOURCE_ROOT="$HOME/nekStab" #--> path to nekStab
#export NEK_SOURCE_ROOT="$NEKSTAB_SOURCE_ROOT/Nek5000"
#export PATH=$NEK_SOURCE_ROOT/bin:$PATH

export USR="matvec.o \
    krylov_decomposition.o \
    eigensolvers.o \
    lapack_wrapper.o \
    fixedp.o \
    usr_extra.o \
    utils.o \
    IO.o \
    postproc.o \
    newton_krylov.o \
    sensitivity.o \
    "

echo "include ${NEKSTAB_SOURCE_ROOT}/core/makefile_nekStab" > makefile_usr.inc
cp --verbose ${NEKSTAB_SOURCE_ROOT}/core/NEKSTAB NEKSTAB.inc
export SOURCE_ROOT=$NEK_SOURCE_ROOT
export N_VERSION=$(cd $NEK_SOURCE_ROOT && git describe --tag --long --always)
echo 'Nek5000 version: '$N_VERSION
export NS_VERSION=$(git describe --tag --long --always)
echo 'nekStab version: '$NS_VERSION

export FFLAGS+=" -DNVERSION=\"'${N_VERSION}'\" -DNSVERSION=\"'${NS_VERSION}'\""

function error_quit {
    echo -e "$@"
    echo
    echo -e 'Usage:'
    echo -e './compile_script --clean'
    echo -e '   To clean build direcrtory. Makenek will ask for cleaning 3rd-party libraries.'
    echo
    echo -e './compile_script --all'
    echo -e '   To compile the code.'
    exit 1
}
