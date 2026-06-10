# Cylinder BoostConv — Base flow via BoostConv acceleration

## Physics
2D flow around a circular cylinder at Re = 100. Computes the unstable steady
base flow using the BoostConv Krylov-accelerated time-stepping method.

## nekStab Mode
`userParam01 = 1.2` — BoostConv

## Prerequisites
- Initial condition: `BF_seed_1cyl0.f00001` (canonical 2128-element seed shared
  with other cylinder baseflow cases)

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

- 6365 BoostConv steps to |F| = 9.93e-10 at t = 335.22
- Wall time 2780 s (≈ 46 min)
- Outputs: converged `BF_1cyl0.f00001` (5.3 MB), `residu.dat`

## Reference
Citro et al. (2017), J. Fluid Mech. 813.
