# Cylinder Animate Modes — Re=100

## Physics
2D flow around a circular cylinder at Re = 100. Reconstructs and animates the
leading direct eigenmode over one oscillation period (T = 2π/ω ≈ 8.56 with
ω = 0.7336851 from `../310_stability_direct/direct/`).

## nekStab Mode
`userParam01 = 4.5` — Animate modes (period-step reconstruction)

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (tracked here; the base flow of `../310_stability_direct/direct/`)
- Direct eigenmodes: `dRe1cyl0.f00002`, `dIm1cyl0.f00002` (tracked here; the eigenmode of `../310_stability_direct/direct/`)
- Eigenvalue file `Spectre_NSd_conv.dat`: not tracked. The routine reads the growth rate and the frequency from it and stops with `ABORT: Could not open Spectre_NSd_conv.dat` when it is missing. Run `../310_stability_direct/direct/` (about 4 min on 8 ranks) and copy `Spectre_NSd_conv.dat` from there.

## Mesh
2128 elements (`lelg=2128` in `SIZE`). Same mesh as the rest of the cylinder
isothermal suite.

## Run
```bash
mks 1cyl
printf "1cyl\n%s/\n" "$(pwd)" > SESSION.NAME   # nek5000 reads the case name from here
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Outputs (regenerated per run)
- `dQ_1cyl0.f00001` … `dQ_1cyl0.f00092` — eigenmode snapshots over one period (spacing 0.0931 = T/92, first snapshot at 0.0931)
- `dQ_1cyl.nek5000` — collection file for ParaView

## Latest Result
Verified 2026-10-09 on 8 ranks:

- Run completes in about 4 s (pure field reconstruction)
- 92 mode snapshots generated covering one full period (last time stamp 8.5639)

## Plot
```bash
python3 plot.py
```

Visualizes the period-animated mode for the canonical Re=100 von Kármán
shedding instability.

## Reference
`ref/reference.json` holds regression values of this repository, not literature values: the time stamp
of the last snapshot (8.5639, one full period), and the maximum of `vx` and the norm of `pr` in the first
snapshot. They follow from the unit-norm direct mode of `../310_stability_direct/direct`. Check a run with
`uv run --with numpy --with pymech python scripts/check_against_ref.py example/cylinder_re100/410_postproc_animate_modes`.
