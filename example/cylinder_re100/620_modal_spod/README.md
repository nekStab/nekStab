# Cylinder Re=100 — SPOD modal analysis (uparam01=6.2)

## What this case shows

Spectral Proper Orthogonal Decomposition yields frequency-resolved spatial
modes that are optimal in the POD sense at each frequency. With `spod_nfft=64`
and `spod_noverlap=48` (Welch's method), the frequency resolution is
df = 1/(64 * 0.5) = 0.03125 Strouhal units and 3 spectral blocks are formed
from 100 snapshots. The leading SPOD mode at Re=100 captures vortex shedding
at St ~= 0.164-0.167.

## Bootstrap IC

`rst_1cyl0.f00001` — copy from `000_dns/rst_1cyl0.f00001` after the DNS
reaches a statistically periodic state (t > 150). This is the restart file
from which the 100 snapshots are collected starting at writeInterval=0.5.

## Expected outputs

- `sRe1cyl0.f00001` ... `sRe1cyl0.f00009` — leading SPOD modes (real part)
- `sIm1cyl0.f00001` ... `sIm1cyl0.f00009` — leading SPOD modes (imag part)
- `spod_spectrum.dat`, `spod_spectrum.png` — energy spectrum vs frequency
- `mea1cyl0.f00001` — time-averaged (mean) flow field
- `logfile` — Slurm job log

## How to run

```bash
mks 1cyl
sbatch run.local.slurm
```

The job reads `1cyl0.f00001` through `1cyl0.f00100` from the working directory
(`modal_prefix='   '`, i.e. SESSION name, 100 snapshots at dt=0.5).
