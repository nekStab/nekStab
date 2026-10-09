# cubic_cavity_re1914 / 000_dns

**uparam01**: 0 (DNS), cold start, Re = 1914, 16 ranks

DNS of the cubic cavity from rest. The flow reaches a steady state by
t = 100 (the probe signal changes by 1e-8 per time unit; 27 min, 69 571 steps on 16 ranks). The last
field, `cav0.f00010`, is the seed of `../210_baseflow_newton/fp/`.

Run: `sbatch run.local.slurm`.

## Reference
`ref/reference.json` checks the last snapshot, `cav0.f00010`: its time stamp (t = 100), the maximum and the norm of `vx` (1 %) and the norm of `vy` (2 %).
