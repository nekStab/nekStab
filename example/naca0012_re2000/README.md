# NACA 0012 — Newton-GMRES base flow

## Physics
2D flow around a NACA 0012 airfoil at Re = 2000. Computes the steady base flow using Newton-GMRES with sponge regions.

## nekStab Mode
`userParam01 = 2.0` — Newton-GMRES for fixed points
- `userParam07 = 220` — Krylov subspace dimension
- Sponge: left = 1, right = 5, strength = 1.7

## Prerequisites
- Initial guess: `BF_naca00120.f00001`

## Run
```bash
mks naca0012        # Compile
nekbmpi naca0012 N  # Run on N MPI ranks
```

## Mesh
4368 elements (`lelg=4368` in `SIZE`, `lpmin=8`).  Polynomial order N=5
(`lx1=6`).

## Run
```bash
mks naca0012
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Latest Result
Verified 2026-05-17 on 8 ranks (Slurm):

- 10 Newton outer iterations to |F| = 9.4e-10
- 27957 inner solver steps, 194 cumulative GMRES iters
- Wall time 722 s (≈ 12 min)
- Outputs: converged `BF_naca00120.f00001` (6.3 MB)

## References

**Methodology (global stability + sensitivity on thin-airfoil wakes):**
Busquet, D., Marquet, O., Richez, F., Juniper, M., & Sipp, D. (2021).
*Global stability, sensitivity and passive control of low-Reynolds-number
flows around NACA 4412 swept wings.* Journal of Fluid Mechanics, 928, A24.
[doi:10.1017/jfm.2021.823](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/global-stability-sensitivity-and-passive-control-of-lowreynoldsnumber-flows-around-naca-4412-swept-wings/9A18EB3560D675EC11D91A7BD4E302F3)

The reference is for the *methodology* — global stability, structural
sensitivity, and passive control on a low-Reynolds-number thin-airfoil
wake. Their case is NACA 4412 swept wings, ours is NACA 0012 (no sweep),
but the matrix-free time-stepper + Krylov-Schur + adjoint + wavemaker
pipeline is the same one nekStab implements here.

## Note on operating points

The active case targets **Re = 2000** only. A `naca0012_Re2500.png` figure
exists in `validation/figures/` as a historical artifact from a higher-Re
run but is not catalog-referenced; the on-disk `.par` is at Re=2000 and
that's the canonical operating point. Re=2500 evidence is treated as
stale until a matching `.par` is re-introduced.
