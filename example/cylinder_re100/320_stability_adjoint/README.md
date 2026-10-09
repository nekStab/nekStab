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
- Base flow: `BF_1cyl0.f00001` (committed: the Newton base flow of
  `../210_baseflow_newton/fp`, copied here; the same file as in `../310_stability_direct/direct`)

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Local I/O Warning
On the earlier 1996-element mesh, too many MPI ranks for local runs corrupted the
`elmap` in generated field files. The runs on the 2128-element mesh use **8 ranks** (`lpmin = 8`).

## Expected Output
Leading adjoint eigenpair σ = 0.1248 ± 0.7337i — **identical** to the direct
spectrum (the operator and its adjoint share spectrum). The spatial mode
differs: concentrated upstream of and around the cylinder rather than
downstream in the wake.

## Mesh
2128 elements (`lelg=2128` in `SIZE`).  Same mesh as the rest of the cylinder
isothermal suite.

## Latest Result
Gate run of 2026-10-09 (Jean Zay, `ifort`, 8 ranks), 2128 mesh, base flow of `../210_baseflow_newton/fp`:

- 2 eigenvalues converged with residual <= 6.1e-7
- Leading mode: sigma = 0.1247881, omega = 0.7336846 (matches the direct sigma to 6e-7)
- Wall time: 229 s
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
