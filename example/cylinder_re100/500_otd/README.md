# Cylinder Re=100 — OTD modes (uparam01=5.0)

## What this case shows

Optimally Time-Dependent (OTD) modes track the instantaneous dominant
instability subspace as the base flow evolves in time. Unlike global
eigenmodes, OTD modes capture transient and time-varying dynamics by
solving a coupled set of linearised equations alongside the nonlinear DNS.
At Re=100 the OTD basis converges to the leading Lyapunov exponents and
their associated spatial structures.

**Note:** `lpert=4` in this scaffold ships a 4-mode OTD basis. Set
`lpert=N` in `SIZE` and recompile with `mks 1cyl` for a larger basis.

## Bootstrap IC

`rst_1cyl0.f00001` — copy from `600_modal_pod/rst_1cyl0.f00001` (DNS
limit-cycle snapshot, lx1=6). OTD evolves alongside the unsteady DNS;
the perturbation basis is initialised with white noise on top of this
DNS field. This matches the canonical OTD usage (time-evolving base
flow), not the variant that linearises around a steady baseflow.

## Expected outputs

- `otd_growth_rates.dat` — instantaneous Lyapunov exponents vs time
- `otd_eigenvalues.dat`  — OTD reduced operator eigenvalues vs time
- `otd_residuals.dat`    — orthonormality residuals
- `dQ_1cyl0.f*****`     — OTD mode field snapshots (every writeInterval=5)
- `logfile`             — Slurm job log

## How to run

```bash
mks 1cyl
sbatch run.local.slurm
```
