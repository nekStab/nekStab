# Cylinder RANS — Reynolds-averaged Navier-Stokes

## Physics
2D flow around a circular cylinder at Re = 40,000. RANS simulation with turbulence model (k-tau), stress formulation, and OIFS time stepping (CFL = 2.5).

## nekStab Mode
`userParam01 = 0` — DNS (RANS via variable properties)
- `variableProperties = yes`, `stressFormulation = yes`
- Scalars: SCALAR01 (tke), SCALAR02 (tau)

## Prerequisites
- Optional restart: `BF_1cyl0.f00001` (commented out; cold start by default)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
RANS solution at Re = 40,000. Turbulent kinetic energy and dissipation rate fields.

## Reference
Internal example.
