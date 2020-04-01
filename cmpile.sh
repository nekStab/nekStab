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
export SOURCE_ROOT="$HOME/Nek5000" #--> path to main code

## INTEL
export FC="mpiifort"; export CC="mpiicc"
export FFLAGS="-mcmodel=medium -shared-intel -xHost -extend-source -mkl -g -traceback"
#export USR_LFLAGS+="-I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl"

## GCC #sudo apt install -y libmpich-dev mpich libblas-dev liblapack-dev
#export FC="mpif77"; export CC="mpicc" #sudo apt install -y libmpich-dev mpich
#export FFLAGS="-mcmodel=medium -march=native -ffixed-line-length-none -g -fbacktrace"
#export USR_LFLAGS+="-L/usr/lib -lblas -llapack" #sudo apt install -y libblas-dev liblapack-dev

export USR="x_eigensolvers.o x_linalg.o x_fixed_point.o x_usr_extra.o x_utilities.o x_IO.o x_postprocessing.o"
#export PPLIST="VENDOR_BLAS" # this substitute the main code eigensolvers

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
