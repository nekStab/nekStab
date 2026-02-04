#!/bin/bash
# Test safe ifx optimization flags (avoiding IPO which breaks with external libs)

set -e

CASE="1cyl"
NPROCS=8
RESULTS_DIR="flag_tests"
NEKSTAB_ROOT="${NEKSTAB_SOURCE_ROOT:-$HOME/nekStab}"
NEK_ROOT="${NEK_SOURCE_ROOT:-$NEKSTAB_ROOT/Nek5000}"

mkdir -p $RESULTS_DIR

# Source Intel environment
source /opt/intel/oneapi/setvars.sh --force 2>/dev/null

clean_build() {
    rm -rf obj nek5000 drive.o makefile makefile_usr.inc NEKSTAB.inc .state .usr 2>/dev/null || true
    rm -rf $NEK_ROOT/3rd_party/gslib/lib 2>/dev/null || true
}

run_test() {
    local name=$1
    local extra_flags=$2

    echo ""
    echo "=============================================="
    echo "Testing: $name"
    echo "Extra flags: $extra_flags"
    echo "=============================================="

    clean_build

    export NEKSTAB_FC=ifx
    export NEKSTAB_EXTRA_FFLAGS="$extra_flags"

    # Build
    echo "Building..."
    if ! mks $CASE > $RESULTS_DIR/${name}_build.log 2>&1; then
        echo "BUILD FAILED for $name"
        tail -10 $RESULTS_DIR/${name}_build.log
        echo "$name,FAILED,FAILED,FAILED,FAILED" >> $RESULTS_DIR/results.csv
        unset NEKSTAB_EXTRA_FFLAGS
        return
    fi

    if [ ! -f nek5000 ]; then
        echo "BUILD FAILED - no nek5000 for $name"
        echo "$name,FAILED,FAILED,FAILED,FAILED" >> $RESULTS_DIR/results.csv
        unset NEKSTAB_EXTRA_FFLAGS
        return
    fi

    # Setup run
    rm -f dRe* dIm* Spectre*.dat logfile 2>/dev/null || true
    echo "$CASE" > SESSION.NAME
    echo "$(pwd)/" >> SESSION.NAME

    # Run
    echo "Running..."
    local start=$(date +%s.%N)
    mpirun -np $NPROCS ./nek5000 > $RESULTS_DIR/${name}_run.log 2>&1
    local wall=$(echo "$(date +%s.%N) - $start" | bc)

    # Extract results
    local nek_time=$(grep "total elapsed time" $RESULTS_DIR/${name}_run.log | awk '{print $(NF-1)}')
    local sigma=$(grep "sigma" $RESULTS_DIR/${name}_run.log | head -1 | awk '{print $3}')
    local omega=$(grep "omega" $RESULTS_DIR/${name}_run.log | head -1 | awk '{print $3}')

    echo "$name,$nek_time,$wall,$sigma,$omega" >> $RESULTS_DIR/results.csv
    echo "  Time: ${nek_time}s, σ=$sigma, ω=$omega"

    unset NEKSTAB_EXTRA_FFLAGS
}

# Initialize results file
echo "name,nek_time,wall_time,sigma,omega" > $RESULTS_DIR/results.csv

echo "Starting flag optimization tests (safe flags only)..."

# Test 1: Baseline (current ifx, no extra flags)
run_test "baseline" ""

# Test 2: -align array64byte (better AVX-512 alignment)
run_test "align64" "-align array64byte"

# Test 3: -funroll-loops (loop unrolling)
run_test "unroll" "-funroll-loops"

# Test 4: -qopt-zmm-usage=high (use more 512-bit registers)
run_test "zmm_high" "-qopt-zmm-usage=high"

# Test 5: -no-prec-div (faster division, slight accuracy loss)
run_test "no_prec_div" "-no-prec-div"

# Test 6: Combination of safe flags
run_test "safe_combo" "-align array64byte -funroll-loops -qopt-zmm-usage=high"

# Test 7: All including no-prec-div
run_test "all_safe" "-align array64byte -funroll-loops -qopt-zmm-usage=high -no-prec-div"

# Test 8: Try -O3 instead of -O2 (more aggressive but safe)
# Note: This requires modifying the build system, skip for now

# Print summary
echo ""
echo "=============================================="
echo "RESULTS SUMMARY"
echo "=============================================="

python3 << 'EOF'
import csv

ref_sigma = 0.1563731e-01  # GCC reference
ref_omega = 0.7565480e+00

print(f"\n{'Config':<15} {'Time (s)':>12} {'Speedup':>10} {'Δσ':>14} {'Δω':>14}")
print("-" * 70)

baseline_time = None
with open('flag_tests/results.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row['name']
        if row['nek_time'] == 'FAILED':
            print(f"{name:<15} {'FAILED':>12}")
            continue

        try:
            nek_time = float(row['nek_time'])
        except:
            print(f"{name:<15} {'PARSE ERROR':>12}")
            continue

        if baseline_time is None:
            baseline_time = nek_time

        speedup = baseline_time / nek_time

        try:
            sigma = float(row['sigma'])
            omega = float(row['omega'])
            d_sigma = sigma - ref_sigma
            d_omega = omega - ref_omega
            print(f"{name:<15} {nek_time:>12.2f} {speedup:>9.2f}x {d_sigma:>+14.2e} {d_omega:>+14.2e}")
        except:
            print(f"{name:<15} {nek_time:>12.2f} {speedup:>9.2f}x {'N/A':>14} {'N/A':>14}")

print()
print("Speedup relative to ifx baseline (no extra flags)")
print("Δσ/Δω relative to GCC reference values")
EOF

echo ""
echo "Results saved in: $RESULTS_DIR/"
