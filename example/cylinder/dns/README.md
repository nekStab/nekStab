# Cylinder DNS — Direct numerical simulation

## Physics
2D flow around a circular cylinder at Re = 50. Vortex shedding develops from a restart file.

## nekStab Mode
`userParam01 = 0` — DNS

## Prerequisites
- Restart file: `rst_1cyl0.f00001`

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
DNS integration to t = 100. Produces field files showing vortex shedding dynamics.

## Reference
Barkley & Henderson (1996), JFM 322.
