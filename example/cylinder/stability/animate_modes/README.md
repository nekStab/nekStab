# Cylinder Animate Modes — Re=100

## Physics
2D flow around a circular cylinder at Re = 100. Reconstructs and animates the
leading direct eigenmode over one oscillation period (T = 2π/ω ≈ 8.56 with
ω = 0.7337 from `../direct/`).

## nekStab Mode
`userParam01 = 4.5` — Animate modes (period-step reconstruction)

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (from `../../baseflow/newton/`)
- Direct eigenmodes: `dRe1cyl0.f0000{1,2}`, `dIm1cyl0.f0000{1,2}` (from `../direct/`)

## Mesh
2128 elements (`lelg=2128` in `SIZE`). Same mesh as the rest of the cylinder
isothermal suite.

## Run
```bash
mks 1cyl
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Outputs (regenerated per run)
- `dQ_1cyl0.f00001` … `dQ_1cyl0.f00010` — eigenmode snapshots over one period
- `dQ_1cyl.nek5000` — collection file for ParaView

## Latest Result
Verified 2026-05-16 on 8 ranks (Slurm):

- Run completes in 0.77 s (pure field reconstruction)
- 10 mode snapshots generated covering one full period

## Plot
```bash
python3 plot.py
```

Visualizes the period-animated mode for the canonical Re=100 von Kármán
shedding instability.

## Reference
Internal example.
