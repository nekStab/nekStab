#!/bin/bash
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

export CASE="1cyl" #--> case name goes here
export NEKSTAB_SOURCE_ROOT="$HOME/nekStab" #--> path to nekStab
export SOURCE_ROOT="$NEKSTAB_SOURCE_ROOT/Nek5000" #--> path to main code
export PATH=$NEKSTAB_SOURCE_ROOT/Nek5000/bin:$PATH

## GCC
#export FC="mpif77"; export CC="mpicc" #sudo apt install -y libmpich-dev mpich
#export FFLAGS="-mcmodel=medium -march=native -ffixed-line-length-none -g -fbacktrace"
#export USR_LFLAGS+="-L/usr/lib -lopenblas" #sudo apt install -y libopenblas-dev

## INTEL
export FC="mpiifort"; export CC="mpiicc"
export FFLAGS="-mcmodel=medium -shared-intel -xHost -extend-source -mkl -g -traceback"
#export USR_LFLAGS+="-I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl" #optional link for old versions

#export PPLIST="VENDOR_BLAS" # this substitute the main code eigensolvers
export USR="x_eigensolvers.o \
	    x_linalg.o \
	    x_fixed_point.o \
	    x_usr_extra.o \
	    x_utilities.o \
	    x_IO.o \
	    x_postprocessing.o"

echo "include ${NEKSTAB_SOURCE_ROOT}/core/makefile_nekStab.inc" > makefile_usr.inc

args=("$@")
argsnr=$#

# check arguments
# parameters number check
if [ $[argsnr] -ne 1 ]; then
    error_quit 'Wrong arguments number!'
fi

for il in "$@"
do
case $il in
        clean)
                ${SOURCE_ROOT}/bin/makenek clean
                shift
                ;;
        all)
                ${SOURCE_ROOT}/bin/makenek ${CASE}
                shift
                ;;
esac
done
