#!/bin/bash
# Benchmark script: Compare GCC vs ifort vs ifx performance and accuracy
# Runs the cylinder DNS test case with all three compilers

set -e

CASE="1cyl"
NPROCS=8          # 8 P-cores (i7-14700: 8P + 12E hybrid)
ENDTIME=100       # Simulation end time
RESULTS_DIR="benchmark_results"

# Pin MPI to P-cores only (CPUs 0-15, avoid E-cores 16-27)
export I_MPI_PIN_DOMAIN=auto
export I_MPI_PIN_PROCESSOR_LIST=0-15
export OMP_NUM_THREADS=1
NEKSTAB_ROOT="${NEKSTAB_SOURCE_ROOT:-$HOME/nekStab}"
NEK_ROOT="${NEK_SOURCE_ROOT:-$NEKSTAB_ROOT/Nek5000}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

log()   { echo -e "${GREEN}[BENCH]${NC} $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
info()  { echo -e "${CYAN}[INFO]${NC} $1"; }

clean_build() {
    log "Cleaning build artifacts..."
    rm -rf obj nek5000 drive.o makefile makefile_usr.inc NEKSTAB.inc *.log.* logfile .state .usr 2>/dev/null || true

    # Clean gslib (required when switching compilers)
    if [ -d "$NEK_ROOT/3rd_party/gslib" ]; then
        cd "$NEK_ROOT/3rd_party/gslib"
        ./install clean 2>/dev/null || true
        cd - > /dev/null
    fi
}

build_gcc() {
    log "Building with GCC..."
    clean_build

    # Unset Intel environment
    unset MKLROOT I_MPI_ROOT ONEAPI_ROOT FC CC CXX 2>/dev/null || true
    export PATH=$(echo $PATH | tr ':' '\n' | grep -v intel | tr '\n' ':' | sed 's/:$//')
    export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v intel | tr '\n' ':' | sed 's/:$//')
    export CFLAGS="-O2 -fPIE"

    mks $CASE 2>&1 | tee $RESULTS_DIR/gcc_build.log
    [ -f nek5000 ] || error "GCC build failed"
    cp nek5000 $RESULTS_DIR/nek5000_gcc
    log "GCC build complete ($(du -h nek5000 | cut -f1))"
}

build_ifort() {
    log "Building with ifort (Intel Classic 2024)..."
    clean_build

    # Source Intel 2024.2 environment (has ifort)
    source /opt/intel/oneapi/compiler/2024.2/env/vars.sh 2>/dev/null || error "Intel 2024.2 compiler not found"
    source /opt/intel/oneapi/mkl/2024.2/env/vars.sh 2>/dev/null || true
    source /opt/intel/oneapi/mpi/2021.13/env/vars.sh 2>/dev/null || error "Intel MPI 2021.13 not found"

    which ifort > /dev/null 2>&1 || error "ifort not found"
    info "Using $(ifort --version 2>&1 | grep -v remark | head -1)"

    mks $CASE 2>&1 | tee $RESULTS_DIR/ifort_build.log
    [ -f nek5000 ] || error "ifort build failed"
    cp nek5000 $RESULTS_DIR/nek5000_ifort
    log "ifort build complete ($(du -h nek5000 | cut -f1))"
}

build_ifx() {
    log "Building with ifx (Intel LLVM 2025)..."
    clean_build

    # Source latest Intel oneAPI
    source /opt/intel/oneapi/setvars.sh --force 2>/dev/null || error "Intel oneAPI not found"

    which ifx > /dev/null 2>&1 || error "ifx not found"
    info "Using $(ifx --version 2>&1 | head -1)"

    # Force ifx via environment variable
    export NEKSTAB_FC=ifx
    mks $CASE 2>&1 | tee $RESULTS_DIR/ifx_build.log
    unset NEKSTAB_FC
    [ -f nek5000 ] || error "ifx build failed"
    cp nek5000 $RESULTS_DIR/nek5000_ifx
    log "ifx build complete ($(du -h nek5000 | cut -f1))"
}

run_benchmark() {
    local compiler=$1
    log "Running $compiler benchmark..."

    # Setup
    cp $RESULTS_DIR/nek5000_${compiler} nek5000
    rm -f lift_drag.dat ${CASE}.log.* logfile fort.* 2>/dev/null || true

    # Ensure history points input file exists
    if [ ! -f ${CASE}.his ]; then
        echo -e "1\n0.0 0.0 0.0" > ${CASE}.his
    fi
    echo "${CASE}" > SESSION.NAME
    echo "$(pwd)/" >> SESSION.NAME

    # Set environment for Intel compilers
    if [ "$compiler" == "ifort" ]; then
        source /opt/intel/oneapi/compiler/2024.2/env/vars.sh 2>/dev/null
        source /opt/intel/oneapi/mpi/2021.13/env/vars.sh 2>/dev/null
    elif [ "$compiler" == "ifx" ]; then
        source /opt/intel/oneapi/setvars.sh --force 2>/dev/null
    fi

    # Run with timing (pin to P-cores 0-15)
    local start=$(date +%s.%N)
    if [ "$compiler" == "gcc" ]; then
        # Use OpenMPI explicitly, bind to cores
        /usr/bin/mpirun -np $NPROCS --bind-to core ./nek5000 > $RESULTS_DIR/${compiler}_run.log 2>&1
    else
        # Intel MPI pinning (environment vars set above)
        mpirun -np $NPROCS ./nek5000 > $RESULTS_DIR/${compiler}_run.log 2>&1
    fi
    local wall=$(echo "$(date +%s.%N) - $start" | bc)

    # Extract Nek5000 internal time (handle scientific notation like 6.56278E+00)
    local nek_time=$(strings $RESULTS_DIR/${compiler}_run.log | grep "total elapsed time" | awk '{printf "%.3f", $NF}')

    # Save history points
    cp ${CASE}.his $RESULTS_DIR/${compiler}.his 2>/dev/null || true

    echo "${compiler},${nek_time},${wall}" >> $RESULTS_DIR/timing.csv
    log "$compiler: Nek=${nek_time}s, Wall=${wall}s"
}

# ============== MAIN ==============
echo ""
echo "========================================"
echo "  nekStab Compiler Benchmark"
echo "  GCC vs ifort (Classic) vs ifx (LLVM)"
echo "========================================"
echo ""

mkdir -p $RESULTS_DIR

# Set simulation time
sed -i "s/endTime = .*/endTime = ${ENDTIME}/" ${CASE}.par
info "Simulation: endTime=${ENDTIME}, NPROCS=${NPROCS}"

# Initialize results
> $RESULTS_DIR/timing.csv

# Build all three
echo ""
log "===== PHASE 1: BUILD ====="
build_gcc
build_ifort
build_ifx

# Run all three
echo ""
log "===== PHASE 2: RUN ====="
run_benchmark "gcc"
run_benchmark "ifort"
run_benchmark "ifx"

# Results
echo ""
log "===== RESULTS ====="
echo ""
printf "%-8s %12s %12s %10s\n" "Compiler" "Nek5000 (s)" "Wall (s)" "Speedup"
echo "--------------------------------------------"
gcc_time=$(grep "^gcc," $RESULTS_DIR/timing.csv | cut -d, -f2)
while IFS=, read -r compiler nek wall; do
    speedup=$(echo "scale=2; $gcc_time / $nek" | bc 2>/dev/null || echo "1.00")
    printf "%-8s %12.3f %12.3f %9.2fx\n" "$compiler" "$nek" "$wall" "$speedup"
done < $RESULTS_DIR/timing.csv
echo ""

log "Results saved in $RESULTS_DIR/"
log "Run 'python3 benchmark_plot.py' for plots"
