#!/bin/bash

#add to .bashrc
#export NEKSTAB_SOURCE_ROOT="$HOME/nekStab" #--> path to nekStab
#export NEK_SOURCE_ROOT="$NEKSTAB_SOURCE_ROOT/Nek5000"
#export PATH=$NEK_SOURCE_ROOT/bin:$PATH

export USR="x_eigensolvers.o \
	    x_linalg.o \
	    x_fixed_point.o \
	    x_usr_extra.o \
	    x_utilities.o \
	    x_IO.o \
	    x_postprocessing.o"

echo "include ${NEKSTAB_SOURCE_ROOT}/core/makefile_nekStab.inc" > makefile_usr.inc
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
