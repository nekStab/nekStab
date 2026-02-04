#!/bin/bash
#
# Cylinder Re=100 Modal Analysis Test
#
# Expected results:
#   Vortex shedding at St ≈ 0.164-0.167
#   With D=1, U=1: f ≈ 0.164-0.167 Hz, T ≈ 6.0 time units
#
# Usage:
#   1. Build:     ./run.sh build
#   2. Run DNS:   ./run.sh dns [nprocs]
#   3. Analyze:   ./run.sh modal [nprocs]
#   4. Plot:      ./run.sh plot
#

set -e
CASE=1cyl
NPROCS=${2:-4}

case "$1" in
    build)
        echo "Building nekStab..."
        ../../bin/makeneks
        ;;

    dns)
        echo "Running DNS to generate snapshots..."
        echo "This will take a while (~100 time units at Re=100)"

        # Set mode to DNS
        sed -i 's/userParam01 = .*/userParam01 = 0/' ${CASE}.par

        # Run
        mpirun -np $NPROCS ./nek5000

        echo ""
        echo "DNS complete. Snapshots saved as ${CASE}0.f*"
        echo "Run './run.sh modal' to perform modal analysis"
        ;;

    modal)
        echo "Running modal analysis (POD/DMD/SPOD)..."

        # Count available snapshots
        NSNAP=$(ls -1 ${CASE}0.f* 2>/dev/null | wc -l)
        echo "Found $NSNAP snapshots"

        if [ "$NSNAP" -lt 10 ]; then
            echo "ERROR: Need at least 10 snapshots. Run DNS first."
            exit 1
        fi

        # Update nsnap in .usr if needed
        # (manual step - user should verify modal_nsnap matches)

        # Set mode to Modal Analysis
        sed -i 's/userParam01 = .*/userParam01 = 6/' ${CASE}.par

        # Run
        mpirun -np $NPROCS ./nek5000

        echo ""
        echo "Modal analysis complete!"
        echo "Output files:"
        echo "  pod_energy.dat   - POD eigenvalue spectrum"
        echo "  dmd_spectrum.dat - DMD eigenvalues (growth, frequency)"
        echo "  spod_spectrum.dat - SPOD spectrum (freq vs eigenvalues)"
        echo "  pod0.f*          - POD modes"
        echo "  dRe0.f*, dIm0.f* - DMD modes (real/imag)"
        echo "  sRe0.f*, sIm0.f* - SPOD modes (real/imag)"
        ;;

    plot)
        echo "Plotting results..."
        python3 plot_modal.py
        ;;

    clean)
        echo "Cleaning..."
        rm -f *.f????? nek5000 *.log logfile *.his
        rm -f pod_energy.dat dmd_spectrum.dat spod_spectrum.dat
        rm -f makefile obj/*.o
        ;;

    *)
        echo "Usage: $0 {build|dns|modal|plot|clean} [nprocs]"
        echo ""
        echo "  build       - Compile nekStab"
        echo "  dns [n]     - Run DNS to generate snapshots (n processors)"
        echo "  modal [n]   - Run modal analysis on snapshots"
        echo "  plot        - Plot results with Python"
        echo "  clean       - Remove output files"
        exit 1
        ;;
esac
