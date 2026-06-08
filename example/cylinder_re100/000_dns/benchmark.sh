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
    rm -rf obj nek5000 drive.o makefile makefile_usr.inc NEKSTAB.inc *.log.* logfile .state .usr 2>/dev/null || true

    # Clean gslib (required when switching compilers)
    if [ -d "$NEK_ROOT/3rd_party/gslib" ]; then
        (cd "$NEK_ROOT/3rd_party/gslib" && ./install clean 2>/dev/null) || true
    fi
}

clean_intel_env() {
    # Remove all Intel paths to avoid MPI version conflicts
    unset MKLROOT I_MPI_ROOT ONEAPI_ROOT FC CC CXX FI_PROVIDER_PATH 2>/dev/null || true
    export PATH=$(echo $PATH | tr ':' '\n' | grep -v intel | tr '\n' ':' | sed 's/:$//')
    export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v intel | tr '\n' ':' | sed 's/:$//')
}

setup_gcc_env() {
    clean_intel_env
    export CFLAGS="-O2 -fPIE"
}

setup_ifort_env() {
    clean_intel_env
    source /opt/intel/oneapi/compiler/2024.2/env/vars.sh 2>/dev/null || error "Intel 2024.2 compiler not found"
    source /opt/intel/oneapi/mkl/2024.2/env/vars.sh 2>/dev/null || true
    source /opt/intel/oneapi/mpi/2021.13/env/vars.sh 2>/dev/null || error "Intel MPI 2021.13 not found"
    command -v ifort > /dev/null 2>&1 || error "ifort not found"
}

setup_ifx_env() {
    clean_intel_env
    source /opt/intel/oneapi/setvars.sh --force 2>/dev/null || error "Intel oneAPI not found"
    command -v ifx > /dev/null 2>&1 || error "ifx not found"
}

build_compiler() {
    local compiler=$1
    local label=$2

    log "Building with $label..."
    clean_build

    case $compiler in
        gcc)   setup_gcc_env ;;
        ifort) setup_ifort_env; info "Using $(ifort --version 2>&1 | grep -v remark | head -1)" ;;
        ifx)   setup_ifx_env; info "Using $(ifx --version 2>&1 | head -1)" ;;
    esac

    [ "$compiler" != "gcc" ] && export NEKSTAB_FC=$compiler
    mks $CASE 2>&1 | tee $RESULTS_DIR/${compiler}_build.log
    unset NEKSTAB_FC

    [ -f nek5000 ] || error "$compiler build failed"
    cp nek5000 $RESULTS_DIR/nek5000_${compiler}
    log "$compiler build complete ($(du -h nek5000 | cut -f1))"
}

run_benchmark() {
    local compiler=$1
    log "Running $compiler benchmark..."

    # Clean up stale processes and sync filesystem
    killall nek5000 2>/dev/null || true
    sync

    # Setup run directory
    cp $RESULTS_DIR/nek5000_${compiler} nek5000
    rm -f lift_drag.dat ${CASE}.log.* logfile fort.* 2>/dev/null || true
    [ -f ${CASE}.his ] || printf "1\n0.0 0.0 0.0\n" > ${CASE}.his
    printf "%s\n%s/\n" "${CASE}" "$(pwd)" > SESSION.NAME

    # Set compiler environment
    case $compiler in
        gcc)   setup_gcc_env ;;
        ifort) setup_ifort_env ;;
        ifx)   setup_ifx_env ;;
    esac

    # Run with timing - verify correct mpirun
    info "Using $(which mpirun) for $compiler"
    local start=$(date +%s.%N)
    if [ "$compiler" = "gcc" ]; then
        /usr/bin/mpirun -np $NPROCS --bind-to core ./nek5000 > $RESULTS_DIR/${compiler}_run.log 2>&1
    else
        mpirun -np $NPROCS ./nek5000 > $RESULTS_DIR/${compiler}_run.log 2>&1
    fi
    local wall=$(echo "$(date +%s.%N) - $start" | bc)

    # Extract Nek5000 internal time (format: "total elapsed time : 5.63728E+01 sec")
    # Use -a to treat log as text (may contain binary from MPI); $5 is the time field
    local nek_time=$(grep -ai "total elapsed time" $RESULTS_DIR/${compiler}_run.log 2>/dev/null | awk '{printf "%.3f", $5}')

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

sed -i "s/endTime = .*/endTime = ${ENDTIME}/" ${CASE}.par
info "Simulation: endTime=${ENDTIME}, NPROCS=${NPROCS}"

> $RESULTS_DIR/timing.csv

echo ""
log "===== PHASE 1: BUILD ====="
build_compiler gcc "GCC"
build_compiler ifort "ifort (Intel Classic 2024)"
build_compiler ifx "ifx (Intel LLVM 2025)"

echo ""
log "===== PHASE 2: RUN (NPROCS=$NPROCS) ====="
run_benchmark gcc
run_benchmark ifort
run_benchmark ifx

echo ""
log "===== RESULTS ====="
echo ""
printf "%-8s %12s %12s %12s\n" "Compiler" "Nek5000 (s)" "Wall (s)" "Speedup"
echo "------------------------------------------------"
slowest_wall=$(sort -t, -k3 -rn $RESULTS_DIR/timing.csv | head -1 | cut -d, -f3)
sort -t, -k3 -rn $RESULTS_DIR/timing.csv | while IFS=, read -r compiler nek wall; do
    speedup_pct=$(echo "scale=4; ($slowest_wall - $wall) / $slowest_wall * 100" | bc 2>/dev/null || echo "0.00")
    printf "%-8s %12.3f %12.3f %+11.2f%%\n" "$compiler" "$nek" "$wall" "$speedup_pct"
done
echo ""

log "Results saved in $RESULTS_DIR/"
log "Run 'python3 benchmark_plot.py' for plots"
