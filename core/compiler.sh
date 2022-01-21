#!/bin/bash

# check presence of Intel compiler in the enviroment
if command -v mpiifort -v 2>/dev/null; then
    echo $(mpiifort -v)
    export FC="mpiifort"
    export CC="mpiicc"
    export FFLAGS+="-extend-source"
    
    if command -v $(ifort --version | grep ^oneapi) 2>/dev/null; then
        export FFLAGS+=""
        #export FFLAGS+=" -qmkl"
    else
        export FFLAGS+=" -xHost"
        #export FFLAGS+=" -mkl"
    fi
    # #export FFLAGS+=" -g -traceback" # tracer 
    # #export FFLAGS+=" -fpe0 -debug extended" # debugger
    # export USR_LFLAGS+="-I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl" #optional link for old versions

else

    echo "gcc version $(gcc --version | grep ^gcc | sed 's/^.* //g')"
    #echo "$(gfortran --version | grep 'gcc version')"
    export FC="mpif90"
    export CC="mpicc"
    export FFLAGS+="-march=native -ffixed-line-length-none -w -Wall"
    #export FFLAGS+=" -g -fbacktrace" # tracer
    #export FFLAGS+=" -Og -Wextra -Warray-bounds -Warray-temporaries"
    #export FFLAGS+=" -Wconversion -Wno-argument-mismatch -fcheck=all"
    #export FFLAGS+=" -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ggdb -fsanitize=leak -fsanitize=address"
fi

os=""; ar=""
case $(uname -s) in
    Darwin) os="Darwin" 

        case $(uname -p) in
            x86_64) ar="x86_64" 
            #echo "Darwin-x86-64 architecture identified"
            ;;
            arm)  ar="arm" 
            echo "Dawin-arm architecture identified"
            export FFLAGS+=" -mcmodel=small"
            ;;
        esac
        ;;

    Linux)  os="Linux"
    echo "Linux-"$(uname -p) "architecture identified"
    export FFLAGS+=" -mcmodel=large" # Restricts code to 2GB and no memory restriction on data

 ;;
esac

# if [[ $(gfortran -dumpversion | cut -f1 -d.) -gt 10 ]]; then
#     export FFLAGS+=" -fallow-argument-mismatch"
# fi

source ${NEKSTAB_SOURCE_ROOT}/core/nekStab.sh
