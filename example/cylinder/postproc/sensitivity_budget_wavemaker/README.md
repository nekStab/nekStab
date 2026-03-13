# Cylinder Sensitivity, Budget & Wavemaker — Post-processing suite

## Physics
2D flow around a circular cylinder at Re = 50. Computes structural sensitivity (wavemaker), perturbation kinetic energy budget, and related diagnostics from pre-computed direct and adjoint modes.

## nekStab Mode
`userParam01 = 4` — Full post-processing (includes 4.1 + 4.2 + 4.3)

## Prerequisites
- Direct and adjoint eigenmode files (from `../../stability/direct/` and `../../stability/adjoint/`)
- No restart file needed (reads eigenmodes automatically)

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
- Structural sensitivity (wavemaker) field
- Perturbation kinetic energy budget terms (production, dissipation, transport)
- Sensitivity to base flow modifications

## Reference
Giannetti & Luchini (2007), JFM 581.
