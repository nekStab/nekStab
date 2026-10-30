# Cylinder Steady Force Sensitivity — Re=100

## Physics
2D flow around a circular cylinder at Re = 100. Computes the **sensitivity of
the leading eigenvalue to a steady body force**, ∇_F λ. The magnitude
|∇_F λ| identifies the regions where a small steady forcing most effectively
modifies the von Kármán instability.

## nekStab Mode
`userParam01 = 4.41` — Sensitivity to steady forcing
- `userParam07 = 50` — Krylov subspace size

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
Run 2026-10-30 on 4 ranks (gfortran, source `dev` at 89770900), T = 1, `k_dim` = 50:

- The GMRES solve (`residu_gmres.dat`) ends at a squared relative residual of
  9.2e-7 after 2 restarts (target 1e-6). The reference holds this as a gate:
  a fresh run must end below 1e-6.
- The sensitivity field `fsr1cyl0.f00001` is concentrated in the recirculation
  bubble and the shear layers next to the cylinder.
- The velocity norms of the field agree to 3e-4 between T = 1 and T = 2, as
  they must: the solution of `(I - exp(T L)) q = r` does not depend on T.
- The `.usr` turns only the outflow `O` into a Dirichlet condition `v`, as in
  stage 320. Turning the inflow `v` into `O` as well made GMRES stall at 1e-2.
- The reference values are regression values of this repository, not literature
  values.

## Plot
```bash
python3 plot.py
```

Shows the velocity magnitude of the sensitivity field (`fsr1cyl0.f00001`).

## Reference
Marquet, Sipp & Jacquin (2008), J. Fluid Mech. 615.
