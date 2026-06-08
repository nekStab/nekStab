# Cylinder DNS — Re=100

## Physics
2D flow around a circular cylinder at Re = 100 (viscosity = -100.0).
Direct numerical simulation; the von Kármán vortex street develops from the
restart file.

## nekStab Mode
`userParam01 = 0` — DNS

## Prerequisites
- Restart file: `rst_1cyl0.f00001` (canonical 2128-element starting field;
  also serves as the universal `BF_seed_1cyl0.f00001` source for the
  cylinder baseflow cases)

## Mesh
2128 elements (`lelg=2128` in `SIZE`, `lpmin=8`). Canonical cylinder mesh
shared across the isothermal cylinder suite.

## Run
```bash
mks 1cyl
sbatch run.local.slurm         # 8 MPI ranks via Slurm
qs
tail -f logfile
```

The bundled `run.local.slurm` writes stdout/stderr to the case-local `logfile`.

Smoke-test the Slurm/MPI placement without running the DNS:

```bash
sbatch --export=ALL,NEK_DRY_RUN=1 run.local.slurm
```

## Latest Result
Verified 2026-05-17 on 8 ranks (Slurm):

- endTime = 100 reached
- Wall time 76 s
- Outputs: `1cyl0.f00001` (final state at t=100), `lift_drag.dat`

The lift coefficient oscillation visible in `lift_drag.dat` shows the
established vortex shedding signature.

## Reference
Barkley & Henderson (1996), J. Fluid Mech. 322.
