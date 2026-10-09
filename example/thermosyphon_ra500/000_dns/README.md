# thermosyphon_ra500 / 000_dns

**Stage**: thermal DNS (userParam01=0), cold start. Buoyancy-driven loop; settings mirror the 310 stage (ifheat, conductivity).

Run: `sbatch run.local.slurm` (4 ranks; 2 min on Jean Zay with `ifort`, 2026-10-09).
`plot.py` writes `plot_snapshot.png` (temperature of the last snapshot); `ref/` holds a copy as `snapshot.png`.

## Reference
`ref/reference.json` checks the last snapshot, `tsyphon0.f00002`: its time stamp (t = 200), the maximum of `vx` (5 %) and the norm of the temperature (2 %). The check is loose.
