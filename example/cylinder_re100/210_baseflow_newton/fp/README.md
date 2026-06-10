# Cylinder Newton — Base flow via Newton-GMRES

## Physics
2D flow around a circular cylinder at Re = 100. Computes the unstable steady
base flow using Newton-GMRES iteration with fixed solver tolerances. Compare
with `newton_dyn/` (Eisenstat-Walker dynamic tolerance scheduling) to see the
reduction in total linear solver calls.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 100` — Krylov subspace dimension (same as newton_dyn)
- `ifdyntol = .false.` — fixed solver tolerances (baseline for comparison)

## Prerequisites
- Initial condition: `BF_seed_1cyl0.f00001` (canonical 2128-element seed shared
  with the rest of the cylinder baseflow suite)

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
Verified 2026-05-15 on 8 ranks (Slurm):

- 5 Newton outer iterations to |F| = 9.83e-10
- 164704 total linear solver steps; 295 cumulative GMRES inner iters
- Wall time 6453 s (≈ 108 min)
- Outputs: converged `BF_1cyl0.f00001` (5.3 MB), `residu_newton.dat`

## Reference
Loiseau et al. (2019), *Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics*.
