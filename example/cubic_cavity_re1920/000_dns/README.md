# cubic_cavity_re1920 / 000_dns

**uparam01**: 0 (DNS), cold start, Re = 1900, 16 ranks

DNS of the cubic cavity at Re = 1900 from rest to t = 100 (21 min, 69 571 steps on 16 ranks,
Jean Zay). Re = 1900 lies below the critical Reynolds number Re_c = 1916.63 of the thesis
(Table 4.1), so the steady flow is stable and the DNS can approach it. At Re = 1920 the steady
flow is unstable and a DNS cannot hold it; the next stage finds it with Newton-GMRES, and it
starts from the last field of this run, `cav0.f00010`.

The flow is close to the steady state but not on it. Between t = 95 and t = 100 the probe signal
still changes by up to 4.6e-4 per time unit (panel (a) of `plot_dns.png`), because the least
damped modes of the steady flow decay slowly. Newton removes the remaining distance.

Run: `sbatch run.local.slurm`.

## Reference
`ref/reference.json` checks the time stamp of the last snapshot `cav0.f00010` (t = 100), the maximum
of `vx`, the norms of `vx` and `vy` (all 1e-3) and the last value of the probe signal `vy` at
(0.25, 0.25, 0) (1e-5). The probe value follows the whole trajectory, including the start-up
oscillation.

## Figure
`python3 plot.py` writes `plot_dns.png`: (a) the probe signal at (0.25, 0.25, 0) and (b) the speed
on the plane z = 0 of the last snapshot. The copy `ref/dns.png` is the figure of the run named in
`ref/reference.json`.
