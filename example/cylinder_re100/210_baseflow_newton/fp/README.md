# Cylinder Newton — Base flow via Newton-GMRES

## Physics
2D flow around a circular cylinder at Re = 100. Computes the unstable steady
base flow using Newton-GMRES iteration with fixed solver tolerances. Compare
with `../dyn/` (Eisenstat-Walker dynamic tolerance scheduling) to see the
reduction in total linear solver calls.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 100` — Krylov subspace dimension (same as `../dyn/`)
- `ifdyntol = .false.` — fixed solver tolerances (baseline for comparison)

## Prerequisites
- Initial condition: `BF_seed_1cyl0.f00001` (canonical 2128-element seed shared
  with the rest of the cylinder baseflow suite)

## Mesh
2128 elements (`lelg=2128` in `SIZE`), polynomial order N=5 (`lx1=6`),
dealiasing N=8 (`lxd=9`).

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Latest Result
Run on 2026-10-09 on 8 ranks of a server (x3, `gfortran`), 5270 s:

- 5 Newton outer iterations to |F| = 9.7e-10 (tolerance 1e-9)
- 93804 time steps in all, 229 cumulative GMRES inner iterations
- Outputs: converged `BF_1cyl0.f00001`, `residu_newton.dat`

The dynamic-tolerance run in `../dyn/` needs 37413 time steps for the same residual.

## Reference
`ref/reference.json` compares the final Newton residual with the tolerance and four statistics of the
converged field (`vx` minimum, maximum and norm, `vy` maximum, tolerance 1e-4). They are the statistics of
the steady flow that `../../110_baseflow_sfd/` reaches by selective frequency damping, an independent
method.

Loiseau et al. (2019), *Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics*.
