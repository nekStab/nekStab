# NACA 0012 — Newton-GMRES base flow

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This airfoil case demonstrates the open-shear-flow matrix-free stability and
wavemaker pipeline on a thin-airfoil wake. The named thesis chapters and
AMR_Krylov examples cite airfoils as part of the broader nekStab application
set, but they do not contain a dedicated NACA 0012 thesis case section. Source
for the method framing: AMR_Krylov_V5/Main `S4_nekstab.tex` and thesis chapter
2, `sec: nekStab`.

## Physics
2D flow around a NACA 0012 airfoil at Re = 2000. Computes the steady base flow using Newton-GMRES with sponge regions.

## nekStab Mode
`userParam01 = 2.0` — Newton-GMRES for fixed points
- `userParam07 = 220` — Krylov subspace dimension
- Sponge: left = 1, right = 5, strength = 1.7

## Prerequisites
- Initial guess of the Newton stage: `BF_seed_naca00120.f00001`, the last snapshot (t = 70) of `000_dns`

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
Run of 2026-10-10 of `210_baseflow_newton/fp`, 8 ranks (Jean Zay, ifort 2021.9):

- Newton from the DNS snapshot: 5 iterations, squared residual 6.08e-7 to 9.49e-10
- 956 s on 8 ranks (220 s on 40 ranks, the same residuals and a base flow that agrees to 9e-12)
- The base flow depends on the seed inside the sponge zone: a second pass seeded with the first output
  moves the field by 3.5e-3, and the growth rate of the direct stage from -0.0643 to -0.0661
  (see `210_baseflow_newton/fp/README.md`)

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

The case runs at **Re = 2000** only: the `.par` file of every stage has
Re = 2000. This release has no Re = 2500 case.
