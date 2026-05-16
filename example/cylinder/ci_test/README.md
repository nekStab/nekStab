# Cylinder CI Test — Minimal integration test

## Physics
2D flow around a circular cylinder at Re = 100. Cold-start DNS used as a fast
build/run sanity check (10 timesteps only).

## nekStab Mode
`userParam01 = 0` — DNS

## Prerequisites
None (cold start, no restart file). Mesh (`1cyl.re2`, `1cyl.ma2`) is a
symbolic link to `../dns/` so this case automatically tracks the canonical
2128-element cylinder mesh.

## Mesh
2128 elements (`lelg=2128` in `SIZE`). Symlinked from `../dns/`.

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Latest Result
Verified 2026-05-17 on 8 ranks (Slurm):

- 10 timesteps, cold-start, no convergence target
- Wall time 0.71 s
- Outputs: `1cyl0.f00001`, `lift_drag.dat`

## Reference
Barkley & Henderson (1996), J. Fluid Mech. 322.
