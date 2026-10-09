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
Gate run of 2026-10-09 (Jean Zay, `ifort`, 8 ranks): the run starts from `rst_1cyl0.f00001` (t = 100) and ends at
`endTime = 200` after 10 593 time steps. Wall time 158 s.

- Outputs: 20 snapshots `1cyl0.f00001` to `1cyl0.f00020` (one every 5 time units), `lift_drag.dat` and the probe
  history `1cyl.his` (v_y at (5, 0)).
- The lift coefficient oscillation in `lift_drag.dat` shows the established vortex shedding.
- `ref/reference.json` checks five quantities. The Strouhal number of the probe signal, 0.1714, is the peak of an FFT
  over the last 70 time units, so its resolution is 1/70 = 0.014 (tolerance 5 %); the DMD stage `../610_modal_dmd`
  gives 0.1657 for the same flow. The time of the last snapshot and the norms of vx and vy and the largest vy in it
  have a tolerance of 1e-6. Jean Zay (`ifort`), `gfortran` at four optimisation levels with and without floating-point
  traps and signalling NaN, and `ifx` at three levels agree to 5e-11 on these norms (`scripts/compiler_battery.py`).
- `plot.py` writes `plot_signal.png`, `plot_fft.png`, `plot_phase.png`, `plot_snapshot.png` and the animation
  `plot.gif`; `ref/` holds copies.

## Reference
Barkley & Henderson (1996), J. Fluid Mech. 322.
