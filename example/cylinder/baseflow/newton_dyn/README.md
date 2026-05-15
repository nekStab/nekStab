# Cylinder Newton Dynamic — Newton-GMRES with dynamic solver tolerances

## Physics
2D flow around a circular cylinder at Re = 100. Computes the unstable steady
base flow using Newton-GMRES with Eisenstat-Walker dynamic solver tolerances
(`ifdyntol = .true.`). Early Newton iterations use relaxed GMRES tolerances
(since the outer residual is large), tightening as Newton converges. This
reduces total linear solver calls relative to the fixed-tolerance Newton case
in `../newton/`.

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
Verified 2026-05-15 on 8 ranks (Slurm):

- 14 Newton outer iterations to |F| = 9.26e-10
- 67357 total linear solver steps, 271 cumulative GMRES inner iters
- Wall time 2510 s (≈ 42 min) -- **~2.6× faster than `../newton/` (108 min)**
- Outputs: converged `BF_1cyl0.f00001` (5.3 MB), `residu_newton.dat`

The dynamic-tolerance speedup is clearer at this mesh resolution: the relaxed
early GMRES iterations skip most of the expensive linear solve work while the
Newton outer residual is still well above 1e-3.

## Reference
Loiseau et al. (2019), *Computational Modelling of Bifurcations and Instabilities in Fluid Dynamics*.
