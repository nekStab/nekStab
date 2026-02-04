#!/bin/bash
# Test different ifx optimization flags for nekStab stability analysis
# Compare performance and eigenvalue accuracy

set -e

CASE="1cyl"
NPROCS=8
RESULTS_DIR="flag_tests"
NEKSTAB_ROOT="${NEKSTAB_SOURCE_ROOT:-$HOME/nekStab}"
NEK_ROOT="${NEK_SOURCE_ROOT:-$NEKSTAB_ROOT/Nek5000}"

# Reference eigenvalues from GCC (our baseline for accuracy)
REF_SIGMA="0.1563731e-01"
REF_OMEGA="0.7565480e+00"

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

    # Export extra flags for makeneks to pick up
    export NEKSTAB_FC=ifx
    export NEKSTAB_EXTRA_FFLAGS="$extra_flags"

    # Build
    echo "Building..."
    mks $CASE > $RESULTS_DIR/${name}_build.log 2>&1

    if [ ! -f nek5000 ]; then
        echo "BUILD FAILED for $name"
        echo "$name,FAILED,FAILED,FAILED,FAILED" >> $RESULTS_DIR/results.csv
        return
    fi

    # Check actual flags used
    grep "FFLAGS:" $RESULTS_DIR/${name}_build.log | head -1 > $RESULTS_DIR/${name}_flags.txt

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

# Test configurations (cumulative - each adds to previous)
echo "Starting flag optimization tests..."

# Baseline: current ifx flags (no fp-model)
run_test "baseline" ""

# Test 1: Add -ipo (interprocedural optimization)
run_test "ipo" "-ipo"

# Test 2: Add -align array64byte (better vectorization alignment)
run_test "align64" "-align array64byte"

# Test 3: Add -funroll-loops
run_test "unroll" "-funroll-loops"

# Test 4: Add -qopt-zmm-usage=high (use more 512-bit registers)
run_test "zmm_high" "-qopt-zmm-usage=high"

# Test 5: Add -no-prec-div (faster division)
run_test "no_prec_div" "-no-prec-div"

# Test 6: Combination of safe flags
run_test "safe_combo" "-ipo -align array64byte -funroll-loops"

# Test 7: All flags combined
run_test "all_flags" "-ipo -align array64byte -funroll-loops -qopt-zmm-usage=high -no-prec-div"

# Print summary
echo ""
echo "=============================================="
echo "RESULTS SUMMARY"
echo "=============================================="
echo ""
echo "Reference (GCC): σ=$REF_SIGMA, ω=$REF_OMEGA"
echo ""

python3 << 'EOF'
import csv

ref_sigma = 0.1563731e-01
ref_omega = 0.7565480e+00

print(f"{'Config':<15} {'Time (s)':>12} {'Speedup':>10} {'Δσ':>14} {'Δω':>14}")
print("-" * 70)

baseline_time = None
with open('flag_tests/results.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row['name']
        if row['nek_time'] == 'FAILED':
            print(f"{name:<15} {'FAILED':>12}")
            continue

        nek_time = float(row['nek_time'])
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
print("Speedup is relative to ifx baseline (no extra flags)")
EOF

echo ""
echo "Detailed results in: $RESULTS_DIR/"
