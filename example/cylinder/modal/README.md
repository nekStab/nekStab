# Cylinder Modal Analysis — POD/DMD/SPOD from DNS snapshots

## Physics
2D flow around a circular cylinder at Re = 100. Modal decomposition of DNS snapshots.

## nekStab Mode
`userParam01 = 6` — Modal analysis (runs POD, DMD, and SPOD together)

Use `6.1`, `6.2`, or `6.3` to run POD, DMD, or SPOD individually.

## Two-Phase Workflow
1. **Phase 1 — Generate snapshots:** set `userParam01 = 0` (DNS). Runs from the restart file and writes field snapshots every `writeInterval = 0.5` time units.
2. **Phase 2 — Modal analysis:** set `userParam01 = 6`. Reads the DNS snapshots and computes POD/DMD/SPOD. The `.par` file is configured for Phase 2.

## Prerequisites
- Restart file: `rst_1cyl0.f00001`
- DNS snapshots from Phase 1 (field files written with `writeInterval = 0.5`)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Modal decomposition results. Expected dominant Strouhal number St ~ 0.164-0.167.

## Reference
Barkley & Henderson (1996), JFM 322.
