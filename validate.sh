#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════════
# nekStab Validation Suite
#
# Runs all 5 AMR paper validation cases (2D only) and compares results against
# expected values with 5% relative error tolerance.
#
# Usage: ./validate.sh [options]
#   --case <name>     Run only specific case (cylinder, backstep, thermosyphon, flipflop, tpjet)
#   --skip-compile    Skip compilation (use existing nek5000 binary)
#   --nprocs <n>      Override number of MPI processes (default: auto-detect P-cores)
#   --dry-run         Show what would be run without executing
#   --help            Show this help message
#
# Ricardo Frantz | Jan 2026
# ═══════════════════════════════════════════════════════════════════════════════

set -e  # Exit on error

# ───────────────────────────────────────────────────────────────────────────────
# Configuration
# ───────────────────────────────────────────────────────────────────────────────
NEKSTAB_ROOT="${NEKSTAB_SOURCE_ROOT:-$HOME/nekStab}"
TOLERANCE=0.05   # 5% relative error
MIN_ELEMS_PER_CORE=20  # Minimum elements per MPI rank for efficiency

# Default to P-cores only on Intel hybrid CPUs
DEFAULT_NPROCS=""
SKIP_COMPILE=false
VALIDATE_ONLY=false
DRY_RUN=false
SINGLE_CASE=""

# ───────────────────────────────────────────────────────────────────────────────
# Colors for terminal output
# ───────────────────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'  # No Color

log_pass() { echo -e "${GREEN}[PASS]${NC} $1"; }
log_fail() { echo -e "${RED}[FAIL]${NC} $1"; }
log_info() { echo -e "${YELLOW}[INFO]${NC} $1"; }
log_detail() { echo -e "${CYAN}      ${NC} $1"; }
log_header() { echo -e "\n${BOLD}═══ $1 ═══${NC}"; }

# ───────────────────────────────────────────────────────────────────────────────
# Results tracking
# ───────────────────────────────────────────────────────────────────────────────
PASSED=0
FAILED=0
declare -a RESULTS=()

# ───────────────────────────────────────────────────────────────────────────────
# Parse command line arguments
# ───────────────────────────────────────────────────────────────────────────────
show_usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --case <name>     Run only specific case"
    echo "                    (cylinder, backstep, thermosyphon, flipflop, tpjet)"
    echo "  --skip-compile    Skip compilation (use existing nek5000 binary)"
    echo "  --validate-only   Only validate existing results (no compile/run)"
    echo "  --nprocs <n>      Override number of MPI processes"
    echo "  --dry-run         Show what would be run without executing"
    echo "  --help            Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                      # Run all cases"
    echo "  $0 --case cylinder      # Run only cylinder case"
    echo "  $0 --validate-only      # Validate existing Spectre_*.dat files"
    echo "  $0 --nprocs 8           # Force 8 MPI processes"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --case)
            SINGLE_CASE="$2"
            shift 2
            ;;
        --skip-compile)
            SKIP_COMPILE=true
            shift
            ;;
        --validate-only)
            VALIDATE_ONLY=true
            shift
            ;;
        --nprocs)
            DEFAULT_NPROCS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help|-h)
            show_usage
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            ;;
    esac
done

# ───────────────────────────────────────────────────────────────────────────────
# Detect P-cores (Performance cores) on Intel hybrid CPUs
# ───────────────────────────────────────────────────────────────────────────────
detect_pcores() {
    local pcores=""

    # Check for Intel hybrid architecture (cpu_core = P-cores, cpu_atom = E-cores)
    if [ -f /sys/devices/cpu_core/cpus ]; then
        pcores=$(cat /sys/devices/cpu_core/cpus)
        echo "$pcores"
    else
        # Non-hybrid CPU: use all cores
        echo "0-$(($(nproc) - 1))"
    fi
}

# Print P-core detection info (call separately from detect_pcores)
print_pcore_info() {
    if [ -f /sys/devices/cpu_core/cpus ]; then
        local pcores=$(cat /sys/devices/cpu_core/cpus)
        log_info "Intel hybrid CPU detected"
        log_detail "P-cores (Performance): $pcores"
        if [ -f /sys/devices/cpu_atom/cpus ]; then
            local ecores=$(cat /sys/devices/cpu_atom/cpus)
            log_detail "E-cores (Efficiency): $ecores (will NOT be used)"
        fi
    else
        log_info "Standard CPU (all cores available)"
    fi
}

# Count CPUs from a range like "0-15" or "0-3,5-7"
count_cpus_in_range() {
    local range="$1"
    local count=0

    # Handle comma-separated ranges
    IFS=',' read -ra RANGES <<< "$range"
    for r in "${RANGES[@]}"; do
        if [[ "$r" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            # Range like "0-15"
            local start="${BASH_REMATCH[1]}"
            local end="${BASH_REMATCH[2]}"
            count=$((count + end - start + 1))
        else
            # Single CPU
            count=$((count + 1))
        fi
    done

    echo "$count"
}

# ───────────────────────────────────────────────────────────────────────────────
# Get optimal core count based on mesh size
# Returns: optimal nprocs (stdout) and info message (stderr)
# ───────────────────────────────────────────────────────────────────────────────
get_optimal_nprocs() {
    local size_file=$1
    local max_cores=$2

    # Extract lelg (total elements) from SIZE file
    local lelg=$(grep -E "parameter\s*\(?\s*lelg\s*=" "$size_file" | sed 's/.*=\s*\([0-9]*\).*/\1/' | head -1)

    if [ -z "$lelg" ]; then
        log_detail "Could not extract lelg from SIZE, using $max_cores cores" >&2
        echo "$max_cores"
        return
    fi

    # Calculate max efficient cores based on mesh size
    local efficient_cores=$((lelg / MIN_ELEMS_PER_CORE))
    [ "$efficient_cores" -lt 1 ] && efficient_cores=1

    # Use minimum of available cores and mesh-limited cores
    if [ "$max_cores" -gt "$efficient_cores" ]; then
        log_detail "Mesh lelg=$lelg too small for $max_cores cores, using $efficient_cores" >&2
        echo "$efficient_cores"
    else
        log_detail "Using $max_cores cores (lelg=$lelg, ${MIN_ELEMS_PER_CORE} elems/core min)" >&2
        echo "$max_cores"
    fi
}

# ───────────────────────────────────────────────────────────────────────────────
# Run a single validation case
# ───────────────────────────────────────────────────────────────────────────────
run_case() {
    local case_name=$1
    local case_dir=$2
    local casename=$3
    local validation_case=$4

    log_header "Case: $case_name"

    local full_dir="$NEKSTAB_ROOT/example/$case_dir"

    if [ ! -d "$full_dir" ]; then
        log_fail "Directory not found: $full_dir"
        FAILED=$((FAILED + 1))
        RESULTS+=("$case_name: FAIL (directory not found)")
        return 1
    fi

    cd "$full_dir"
    log_detail "Working directory: $full_dir"

    # Determine number of processes
    local nprocs
    if [ -n "$DEFAULT_NPROCS" ]; then
        nprocs="$DEFAULT_NPROCS"
        log_detail "Using user-specified $nprocs cores"
    else
        nprocs=$(get_optimal_nprocs "SIZE" "$P_CORE_COUNT")
    fi

    # Get P-core CPU list for taskset binding
    local cpu_list="$P_CORE_CPUS"

    if [ "$DRY_RUN" = true ]; then
        log_info "[DRY RUN] Would compile: mks $casename"
        log_info "[DRY RUN] Would run on $nprocs P-cores (CPUs: $cpu_list)"
        log_detail "Command: mpirun -np $nprocs --bind-to core --cpu-set $cpu_list ./nek5000"
        log_info "[DRY RUN] Would validate: $validation_case"
        RESULTS+=("$case_name: DRY RUN")
        return 0
    fi

    # Validate-only mode: skip compile and run, just validate existing results
    if [ "$VALIDATE_ONLY" = true ]; then
        log_detail "Validate-only mode: checking existing results"
        # Check that output files exist
        local case_info=$(python3 -c "import json; d=json.load(open('$NEKSTAB_ROOT/validation/expected_values.json')); print(d['$validation_case']['output_file'])" 2>/dev/null)
        if [ -z "$case_info" ]; then
            case_info="Spectre_NSd.dat"
        fi
        if [ ! -f "$case_info" ]; then
            log_fail "No results found: $case_info"
            log_detail "Run without --validate-only to generate results"
            FAILED=$((FAILED + 1))
            RESULTS+=("$case_name: FAIL (no results)")
            return 1
        fi
        log_info "Validating existing results..."
        if ! python3 "$NEKSTAB_ROOT/validate_compare.py" "$validation_case" "$full_dir"; then
            log_fail "Validation failed for $case_name"
            FAILED=$((FAILED + 1))
            RESULTS+=("$case_name: FAIL (validation)")
            return 1
        fi
        PASSED=$((PASSED + 1))
        RESULTS+=("$case_name: PASS")
        return 0
    fi

    # Step 1: Compile
    if [ "$SKIP_COMPILE" = false ]; then
        log_info "Compiling $casename..."
        # Ensure nekStab environment is set
        if [ -z "$NEKSTAB_SOURCE_ROOT" ]; then
            export NEKSTAB_SOURCE_ROOT="$NEKSTAB_ROOT"
            export NEK_SOURCE_ROOT="$NEKSTAB_ROOT/Nek5000"
            export PATH="$NEK_SOURCE_ROOT/bin:$NEKSTAB_SOURCE_ROOT/bin:$PATH"
        fi
        if ! mks "$casename" > build.log 2>&1; then
            log_fail "Compilation failed. Check $full_dir/build.log"
            FAILED=$((FAILED + 1))
            RESULTS+=("$case_name: FAIL (compilation)")
            return 1
        fi
        log_detail "Compilation successful"
    else
        log_detail "Skipping compilation (--skip-compile)"
        if [ ! -x "./nek5000" ]; then
            log_fail "No nek5000 binary found and --skip-compile specified"
            FAILED=$((FAILED + 1))
            RESULTS+=("$case_name: FAIL (no binary)")
            return 1
        fi
    fi

    # Step 2: Run simulation with P-core binding
    log_info "Running on $nprocs P-cores (CPUs: $cpu_list)..."

    # Set up SESSION.NAME
    echo "$casename" > SESSION.NAME
    echo "$(pwd)/" >> SESSION.NAME

    # Clean up old files
    rm -f logfile ioinfo Spectre_*.dat

    # Run with CPU pinning to P-cores only
    local run_cmd="mpirun -np $nprocs"

    # Add CPU binding if hwloc/OpenMPI supports it
    if mpirun --help 2>&1 | grep -q "cpu-set"; then
        run_cmd="$run_cmd --bind-to core --cpu-set $cpu_list"
    elif command -v taskset >/dev/null 2>&1; then
        # Fallback to taskset
        run_cmd="taskset -c $cpu_list $run_cmd"
    fi

    run_cmd="$run_cmd ./nek5000"

    log_detail "Command: $run_cmd"

    local start_time=$(date +%s)
    if ! $run_cmd > "${casename}.log.${nprocs}" 2>&1; then
        log_fail "Simulation failed. Check ${casename}.log.${nprocs}"
        FAILED=$((FAILED + 1))
        RESULTS+=("$case_name: FAIL (simulation)")
        return 1
    fi
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))
    log_detail "Simulation completed in ${elapsed}s"

    # Step 3: Validate results
    log_info "Validating results..."

    if ! python3 "$NEKSTAB_ROOT/validate_compare.py" "$validation_case" "$full_dir"; then
        log_fail "Validation failed for $case_name"
        FAILED=$((FAILED + 1))
        RESULTS+=("$case_name: FAIL (validation)")
        return 1
    fi

    PASSED=$((PASSED + 1))
    RESULTS+=("$case_name: PASS")
    return 0
}

# ───────────────────────────────────────────────────────────────────────────────
# Main execution
# ───────────────────────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════════════════════════"
echo "                     nekStab Validation Suite                               "
echo "═══════════════════════════════════════════════════════════════════════════"
echo ""

# Detect P-cores
P_CORE_CPUS=$(detect_pcores)
P_CORE_COUNT=$(count_cpus_in_range "$P_CORE_CPUS")

# Print detection info
print_pcore_info
log_info "Maximum available P-cores: $P_CORE_COUNT (CPUs: $P_CORE_CPUS)"
log_info "Tolerance: ${TOLERANCE} (5%)"
log_info "nekStab root: $NEKSTAB_ROOT"

if [ -n "$DEFAULT_NPROCS" ]; then
    log_info "User-specified cores: $DEFAULT_NPROCS"
fi

if [ "$DRY_RUN" = true ]; then
    log_info "DRY RUN MODE - no commands will be executed"
fi

if [ "$VALIDATE_ONLY" = true ]; then
    log_info "VALIDATE-ONLY MODE - checking existing results only"
fi

echo ""

# Define test cases
# Format: case_name|case_dir|casename|validation_case_name
CASES=(
    "Cylinder 2D (Hopf)|cylinder/stability/direct|1cyl|cylinder_hopf"
    "Backward Step (Transient Growth)|back_fstep/transient_growth|bfs|backstep_tg"
    "Thermosyphon (Pitchfork)|thersyphon/stability/direct|tsyphon|thermosyphon_pitchfork"
    "Flip-flop (Neimark-Sacker)|flip_flop/stability/direct_Floquet|2cyl|flipflop_ns"
    "Forced Jet (Period-doubling)|tpjet/stability/direct_Floquet|tpjet|tpjet_pd"
)

# Filter to single case if specified
if [ -n "$SINGLE_CASE" ]; then
    FILTERED_CASES=()
    for case_entry in "${CASES[@]}"; do
        IFS='|' read -r name dir casename valcase <<< "$case_entry"
        # Match by partial name or validation case name
        if [[ "${name,,}" == *"${SINGLE_CASE,,}"* ]] || [[ "$valcase" == *"$SINGLE_CASE"* ]]; then
            FILTERED_CASES+=("$case_entry")
        fi
    done

    if [ ${#FILTERED_CASES[@]} -eq 0 ]; then
        log_fail "No case matching '$SINGLE_CASE' found"
        echo "Available cases:"
        for case_entry in "${CASES[@]}"; do
            IFS='|' read -r name dir casename valcase <<< "$case_entry"
            echo "  - $name ($valcase)"
        done
        exit 1
    fi

    CASES=("${FILTERED_CASES[@]}")
    log_info "Running only: ${CASES[0]%%|*}"
fi

# Run all cases
TOTAL=${#CASES[@]}
CURRENT=0

for case_entry in "${CASES[@]}"; do
    CURRENT=$((CURRENT + 1))
    IFS='|' read -r name dir casename valcase <<< "$case_entry"

    log_info "Case $CURRENT/$TOTAL: $name"

    run_case "$name" "$dir" "$casename" "$valcase" || true
done

# ───────────────────────────────────────────────────────────────────────────────
# Summary
# ───────────────────────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════════════════════════"
echo "                          VALIDATION SUMMARY                                "
echo "═══════════════════════════════════════════════════════════════════════════"
echo ""

for result in "${RESULTS[@]}"; do
    if [[ "$result" == *"PASS"* ]]; then
        log_pass "$result"
    elif [[ "$result" == *"DRY RUN"* ]]; then
        log_info "$result"
    else
        log_fail "$result"
    fi
done

echo ""
if [ "$DRY_RUN" = true ]; then
    echo -e "${YELLOW}DRY RUN COMPLETE${NC}"
elif [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}${BOLD}ALL $PASSED/$TOTAL CASES PASSED${NC}"
    echo "Ready to merge to main branch."
    exit 0
else
    echo -e "${RED}${BOLD}$FAILED/$TOTAL CASES FAILED${NC}"
    echo "Please investigate failures before merging."
    exit 1
fi
