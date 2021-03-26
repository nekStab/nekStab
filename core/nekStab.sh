#!/bin/bash

#add to .bashrc
#export NEKSTAB_SOURCE_ROOT="$HOME/nekStab" #--> path to nekStab
#export NEK_SOURCE_ROOT="$NEKSTAB_SOURCE_ROOT/Nek5000"
#export PATH=$NEK_SOURCE_ROOT/bin:$PATH

export USR="eigensolvers.o \
    linalg.o \
    fixedp.o \
    usr_extra.o \
    utils.o \
    IO.o \
    postproc.o\
    newton_krylov.o"

echo "include ${NEKSTAB_SOURCE_ROOT}/core/makefile_nekStab" > makefile_usr.inc
cp --verbose ${NEKSTAB_SOURCE_ROOT}/core/NEKSTAB NEKSTAB.inc
export SOURCE_ROOT=$NEK_SOURCE_ROOT

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
