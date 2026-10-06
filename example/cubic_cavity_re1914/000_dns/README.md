# cubic_cavity_re1914 / 000_dns

**uparam01**: 0 (DNS), cold start, Re = 1914, 16 ranks

DNS of the cubic cavity from rest. The flow reaches a steady state by
t = 100 (probe signal changes by 4e-8 per step; 14 min on 16 ranks). The last
field, `cav0.f00010`, is the seed of `../210_baseflow_newton/fp/`.

Run: `sbatch run.local.slurm`.
