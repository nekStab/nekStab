# Cylinder Newton Dynamic — Newton-GMRES with dynamic solver tolerances

## Physics
2D flow around a circular cylinder at Re = 100. Computes the unstable steady
base flow using Newton-GMRES with Eisenstat-Walker dynamic solver tolerances
(`ifdyntol = .true.`). Early Newton iterations use relaxed GMRES tolerances
(since the outer residual is large), tightening as Newton converges. This
reduces total linear solver calls relative to the fixed-tolerance Newton case
in `../fp/`.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 100` — Krylov subspace dimension
- `ifdyntol = .true.` — Eisenstat-Walker dynamic tolerance scheduling

## Prerequisites
- Initial condition: `BF_seed_1cyl0.f00001` (canonical 2128-element seed)

## Mesh
2128 elements (`lelg=2128` in `SIZE`), polynomial order N=7 (`lx1=8`),
dealiasing N=11 (`lxd=12`).

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Latest Result
Run on 2026-10-09 on 8 ranks of Jean Zay (`ifort`), 672 s:

- 14 Newton outer iterations to |F| = 9.3e-10 (tolerance 1e-9)
- 37413 time steps in all, 307 cumulative GMRES inner iterations
- Outputs: converged `BF_1cyl0.f00001`, `residu_newton.dat`

The fixed-tolerance run in `../fp/` needs 93804 time steps for the same residual.

## Reference
`ref/reference.json` compares the final Newton residual with the tolerance and four statistics of the
converged field (`vx` minimum, maximum and norm, `vy` maximum, tolerance 1e-4). They are the statistics of
the steady flow that `../../110_baseflow_sfd/` reaches by selective frequency damping, an independent
method.

Loiseau et al. (2019), *Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics*.
