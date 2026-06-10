# Cylinder Re=100 — DMD modal analysis (uparam01=6.1)

## What this case shows

Dynamic Mode Decomposition extracts spatially coherent modes and their
associated frequencies and growth rates from a snapshot ensemble. Each DMD
mode has a single complex frequency; the leading mode at Re=100 captures the
vortex shedding at St ~= 0.164-0.167. Using `dmd_rank=20` avoids auto-rank
truncation that would miss the subdominant oscillatory modes.

## Bootstrap IC

`rst_1cyl0.f00001` — copy from `000_dns/rst_1cyl0.f00001` after the DNS
reaches a statistically periodic state (t > 150). This is the restart file
from which the 100 snapshots are collected starting at writeInterval=0.5.

## Expected outputs

- `dm11cyl0.f00001`, `dm21cyl0.f00001` — leading DMD mode fields
- `dmd_spectrum.dat` — DMD eigenvalues (frequency, growth rate, amplitude)
- `dmd_svd.dat` — singular values of the snapshot matrix
- `mea1cyl0.f00001` — time-averaged (mean) flow field
- `logfile` — Slurm job log

## How to run

```bash
mks 1cyl
sbatch run.local.slurm
```

The job reads `1cyl0.f00001` through `1cyl0.f00100` from the working directory
(`modal_prefix='   '`, i.e. SESSION name, 100 snapshots at dt=0.5).
