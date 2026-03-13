# Cylinder Newton with Active Scalar — Newton-GMRES with buoyancy coupling

## Physics
2D flow around a circular cylinder at Re = 50 with active temperature field. Buoyancy (`userParam06 = 0.1`) couples temperature to the momentum equation via a vertical body force, deflecting the wake. Newton-GMRES with dynamic solver tolerances (`ifdyntol`) finds the coupled steady base flow.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam06 = 0.1` — buoyancy coefficient (Boussinesq forcing in y)
- `userParam07 = 100` — Krylov subspace dimension (GMRES)
- `ifdyntol = .true.` — dynamic solver tolerances (set in `.usr`)
- Temperature + 3 passive scalars enabled

## Prerequisites
- Optional restart: `BFt_1cyl0.f00001` (commented out; cold start by default)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Converged steady base flow with temperature field. The buoyancy force breaks the top-bottom symmetry of the wake. Newton residual decreases to machine precision.

## Reference
Internal example.
