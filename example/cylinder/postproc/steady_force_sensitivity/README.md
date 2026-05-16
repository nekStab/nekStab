# Cylinder Steady Force Sensitivity — Re=100

## Physics
2D flow around a circular cylinder at Re = 100. Computes the **sensitivity of
the leading eigenvalue to a steady body force**, ∇_F λ. The magnitude
|∇_F λ| identifies the regions where a small steady forcing most effectively
modifies the von Kármán instability.

## nekStab Mode
`userParam01 = 4.41` — Sensitivity to steady forcing
- `userParam07 = 50` — Krylov subspace size (single-node-viable)

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (from `../../baseflow/newton/`)
- Direct modes: `dRe1cyl0.f0000{1,2}`, `dIm1cyl0.f0000{1,2}` (from `../../stability/direct/`)
- Adjoint modes: `aRe1cyl0.f0000{1,2}`, `aIm1cyl0.f0000{1,2}` (from `../../stability/adjoint/`)
- Wavemaker sensitivity inputs: `sr_1cyl0.f00001`, `si_1cyl0.f00001`
  (from `../sensitivity_budget_wavemaker/`)

## Mesh
2128 elements (`lelg=2128` in `SIZE`). Same mesh as the rest of the cylinder
isothermal suite.

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Latest Result
Verified 2026-05-15 on 8 ranks (Slurm):

- Runs an internal adjoint Arnoldi iteration (`residu_arnoldi.dat`) plus a
  GMRES solve (`residu_gmres.dat`) to compute the gradient.
- Wall time 1716 s (≈ 29 min)
- Outputs: `fsr1cyl0.f00001` (steady-force sensitivity field) plus the
  intermediate `sr_/si_` files used as inputs.

## Plot
```bash
python3 plot.py
```

Shows |∇_F λ| modulus identifying the optimal force-placement region.

## Reference
Marquet, Sipp & Jacquin (2008), J. Fluid Mech. 615.
