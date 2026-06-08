# Cylinder Adjoint Stability — Re=100

## Physics
2D flow around a circular cylinder at Re=100. Computes the leading **adjoint**
eigenvalues and modes of the linearized Navier-Stokes operator around the
steady base flow. Adjoint modes peak upstream of the cylinder — the
receptivity region where the flow is most sensitive to perturbations.

## nekStab Mode
`userParam01 = 3.2` — Adjoint eigensolver (Schur-Krylov)
- `userParam07 = 50` — `k_dim`: Krylov subspace size (single-node-viable)
- Sponge: `userParam08 = 5` (left), `userParam09 = 5` (right)

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
`elmap` in generated field files. **8 ranks** is the known-good local count.

## Expected Output
Leading adjoint eigenpair near σ = 0.13 ± 0.74i — **identical** to the direct
spectrum (the operator and its adjoint share spectrum). The spatial mode
differs: concentrated upstream of and around the cylinder rather than
downstream in the wake.

## Mesh
2128 elements (`lelg=2128` in `SIZE`).  Same mesh as the rest of the cylinder
isothermal suite.

## Latest Result
Verified 2026-05-15 on 8 ranks (Slurm), 2128 mesh, BF from `../../baseflow/newton/`:

- 2 eigenvalues converged with residual ≤ 3.0e-7
- Leading mode: σ = 0.1247587, ω = 0.7337421 (matches direct to 6 digits)
- Wall time: 812 s
- Outputs: `aRe1cyl0.f0000{1,2}`, `aIm1cyl0.f0000{1,2}`,
  `Spectre_NSa*.dat`, `Spectre_Ha*.dat`

## Plot

```bash
python3 plot.py
```

Produces a 2-panel figure:
- (a) σ vs f adjoint spectrum from `Spectre_NSa.dat`
- (b) Leading adjoint eigenmode (real part) — note concentration upstream
  of the cylinder

The product of the leading direct and adjoint modes gives the **structural
sensitivity** (wavemaker) field used in the postproc cases.

## Reference
Giannetti & Luchini (2007), J. Fluid Mech. 581.
