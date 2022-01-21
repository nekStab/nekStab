#!/bin/bash

## GCC

# sudo apt install -y libmpich-dev mpich libopenblas-dev

export FC="mpif77"
export CC="mpicc"
export FFLAGS="-march=native -ffixed-line-length-none -w -Wall"

os=""
case $(uname -s) in
    Darwin) os="Darwin" ;;
    Linux)  os="Linux" ;;
esac

ar=""
case $(uname -p) in
    x86_64) ar="x86_64" ;;
    arm)  ar="arm" ;;
esac


if [[ $os == "Linux"* ]]; then
    echo "Compiling in Linux"
    export FFLAGS+=" -mcmodel=large"
    export USR_LFLAGS+="-L/usr/lib" # -lopenblas"

elif [[ "$os" == "Darwin"* ]]; then

    if [[ $ar == "x86_64"* ]]; then
    echo "Compiling for Mac OS-x86-64 architecture"

    elif [[ "$ar" == "arm"* ]]; then
    echo "Compiling for Mac OS-ARM architecture"
    export FFLAGS+=" -mcmodel=small"

fi
fi

if [[ $(gfortran -dumpversion | cut -f1 -d.) -gt 10 ]]; then
    export FFLAGS+=" -fallow-argument-mismatch"
fi


#export FFLAGS+="-g -fbacktrace" # tracer
#export FFLAGS="-Og -Wextra -Warray-bounds -Warray-temporaries -Wconversion -Wno-argument-mismatch -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ggdb -fsanitize=leak -fsanitize=address"


source ${NEKSTAB_SOURCE_ROOT}/core/nekStab.sh
