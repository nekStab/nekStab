# Cylinder CI Test — Minimal integration test

## Physics
2D flow around a circular cylinder at Re = 50. Cold-start DNS used as a quick build/run sanity check.

## nekStab Mode
`userParam01 = 0` — DNS

## Prerequisites
None (cold start, no restart file).

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Runs 10 time steps and exits. Verifies that the code compiles and executes without errors.

## Reference
Barkley & Henderson (1996), JFM 322.
