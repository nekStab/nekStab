# Cylinder Direct Stability — Re=100

## Physics
2D flow around a circular cylinder at Re=100. Computes the leading eigenvalues
of the linearized Navier-Stokes operator around the steady base flow (the
unstable von Kármán shedding instability).

## nekStab Mode
`userParam01 = 3.1` — Direct eigensolver (Schur-Krylov)
- `userParam03 = 2`   — `schur_tgt`: target eigenvalues to converge
- `userParam06 = 0`   — isothermal (no buoyancy coupling)
- `userParam07 = 50`  — `k_dim`: Krylov subspace size (single-node-viable)

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (committed; sourced from
  `../../baseflow/newton/BF_1cyl0.f00001`, the canonical Re=100 Newton BF)

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Local I/O Warning
On this 1996-element mesh, too many MPI ranks for local runs can corrupt the
`elmap` in generated field files (pymech warns "elmap appears to be
corrupted"). **8 ranks** is the known-good local count.

## Expected Output
Leading eigenpair near σ = 0.13 ± 0.74i (growth rate ± angular frequency),
giving Strouhal St = ω/(2π) ≈ 0.118 for the linear instability mode around the
Re=100 unstable steady base flow.  This is the textbook 2D cylinder shedding
mode (Barkley & Henderson 1996).

## Mesh
2128 elements (`lelg=2128` in `SIZE`).  Same mesh as the rest of the cylinder
isothermal suite.

## Latest Result
Verified 2026-05-15 on 8 ranks (Slurm), 2128 mesh, BF from `../../baseflow/newton/`:

- 2 eigenvalues converged with residual ≤ 7.1e-7
- Leading mode: σ = 0.1247590, ω = 0.7337421  (St = 0.1168)
- Wall time: 773 s
- Outputs: `dRe1cyl0.f0000{1,2}`, `dIm1cyl0.f0000{1,2}`, `Spectre_NSd*.dat`,
  `Spectre_Hd*.dat`

The eigenvalue is slightly different from earlier 1996-mesh runs
(σ ≈ 0.1275 + 0.7447i) — the refined 2128 mesh gives a more accurate
linear-stability result.

## Plot

```bash
python3 plot.py
```

Produces a 2-panel figure:
- (a) σ vs f spectrum from `Spectre_NSd.dat`
- (b) Leading direct eigenmode (real part) overlaid on the cylinder geometry

## Reference
Barkley & Henderson (1996), J. Fluid Mech. 322.
