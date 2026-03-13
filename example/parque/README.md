# Wind Farm (Parque) — RANS wind farm simulation

## Physics
Wind farm simulation at Re = 22,500 using RANS with actuator disk model. Includes atmospheric boundary layer profile, surface roughness, and multiple turbine rows with varying thrust coefficients.

## nekStab Mode
DNS (RANS via variable properties)
- `variableProperties = yes`, `stressFormulation = yes`
- Scalars: SCALAR02 (tke), SCALAR03 (tau)

## Prerequisites
None (cold start).

## Run
```bash
mks parque        # Compile
nekbmpi parque N  # Run on N MPI ranks
```

## Expected Output
Steady-state RANS solution of the wind farm wake. Flow fields showing turbine wake interactions.

## Reference
Internal example.
