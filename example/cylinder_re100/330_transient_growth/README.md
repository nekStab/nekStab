# Cylinder Re=100 — transient growth (uparam01=3.3)

## What this case shows

Transient growth computes the optimal perturbation that maximises energy
amplification over a target horizon T around the Newton baseflow. The
Krylov-Schur driver runs paired direct/adjoint time-integrations; internally
uparam(1) is rewritten from 3.3 to 3.1 and back each iteration. At Re=100
the flow is linearly stable, so transient non-normal amplification reveals
mechanisms active before asymptotic decay.

## Bootstrap IC

`BF_1cyl0.f00001` — copy from `210_baseflow_newton/fp/BF_1cyl0.f00001` after
the Newton solver converges (residual < 1e-9).

## Expected outputs

- `Spectre_NSd_conv.dat` — Krylov convergence history
- `Spectre_Hd.dat`, `Spectre_NSd.dat` — Ritz values at each restart
- `dRe1cyl0.f00001`, `dIm1cyl0.f00001` — leading optimal perturbation (real/imag)
- `logfile` — Slurm job log

## How to run

```bash
mks 1cyl
sbatch run.local.slurm
```

Build links the mesh via `1cyl.ma2` / `1cyl.re2` (symlink or copy from
`210_baseflow_newton/fp/`).
