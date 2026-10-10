# cubic_cavity_re1914 / 000_dns

**uparam01**: 0 (DNS), cold start, Re = 1914, 16 ranks

DNS of the cubic cavity from rest to t = 100 (25 min, 69 571 steps on 16 ranks, Jean Zay).
The flow is close to the steady state but not on it. Between t = 95 and t = 100 the probe signal
still changes by up to 5e-4 per time unit (panel (a) of `plot_dns.png`). The least damped modes of
the steady flow have sigma = -1.8e-4 (the T1 pair) and -0.026, so the start-up oscillation decays
slowly and is still present at t = 100. The last field, `cav0.f00010`, is the seed of
`../210_baseflow_newton/fp/`. Newton removes the remaining distance (first residual 9.0e-6).

Run: `sbatch run.local.slurm`.

## Reference
`ref/reference.json` checks the time stamp of the last snapshot `cav0.f00010` (t = 100), the maximum of `vx`, the norms of `vx` and `vy` (all 1e-3) and the last value of the probe signal `vy` at (0.25, 0.25, 0) (1e-5). The probe value follows the whole trajectory, including the start-up oscillation.

## Figure
`python3 plot.py` writes `plot_dns.png`: (a) the probe signal at (0.25, 0.25, 0) and (b) the speed on the plane z = 0 of the last snapshot. The copy `ref/dns.png` is the figure of the run named in `ref/reference.json`.
