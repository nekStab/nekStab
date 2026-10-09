# lid_driven_re3600 / 000_dns — DNS of the regularized lid-driven cavity to the steady flow

**Stage**: 000_dns · **uparam01**: 0 · **Re**: 3600 · **Start from**: cold start (the lid profile in
`useric`) · **Ranks**: 4. Stage rationale: repo-root `EXAMPLES.md`.

## What it does
2D cavity of aspect ratio 1.5 (`userParam10 = 1.5`) with the regularized lid u = (1 − (2x)²)²
(Shen 1991), 100 elements of order 11 (`lx1 = 12`, `lxd = 18`), `filtering = explicit`. Re = 3600
is below the first Hopf bifurcation, so the flow is steady and the DNS relaxes to it. The probe
at (0.25, 0.25) approaches its final value with the rate of the leading eigenvalue, 0.0168, so the
run lasts 600 time units (`endTime`) and writes a snapshot every 100 (`writeInterval`). The last
snapshot, `cav0.f00006`, seeds `../210_baseflow_newton/fp/` (copy it to `BF_cav0.f00001`): at
`lx1 = 12` its first Newton residual is 4.4e-9.

Run: `sbatch run.local.slurm` (4 ranks; 18 min on Jean Zay with `ifort`, 2026-10-09).
`plot.py` writes `plot_snapshot.png` (velocity magnitude of the last snapshot); `ref/` holds a copy as `snapshot.png`.

## Resolution
Order 5 (`lx1 = 6`, the former setting) is too low for this flow: the steady state of its
filtered DNS is not a fixed point of the unfiltered map of the later stages, and they found a
different, unstable state. The resolution table is in `../README.md`.
