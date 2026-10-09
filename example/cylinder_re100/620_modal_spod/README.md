# Cylinder Re=100 — SPOD modal analysis (uparam01=6.3)

## What this case shows

Spectral Proper Orthogonal Decomposition yields frequency-resolved spatial
modes that are optimal in the POD sense at each frequency. With `spod_nfft=64`
and `spod_noverlap=48` (Welch's method), the frequency resolution is
df = 1/(64 * 0.5) = 0.03125 Strouhal units and 3 spectral blocks are formed
from 100 snapshots. The leading SPOD mode at Re=100 captures vortex shedding
at St ~= 0.164-0.167.

## Bootstrap IC

`run.local.slurm` first runs a DNS (`dns.par`, `userParam01 = 0`) from the restart file
`rst_1cyl0.f00001` (a limit-cycle state at t = 150) to t = 200. It writes the 100 snapshots
`1cyl0.f00001` ... `1cyl0.f00100` (one every 0.5 time units), then restores `1cyl.par` and runs the
modal analysis on them. The restart file and `dns.par` are kept in `../600_modal_pod`; the script
copies them when they are missing. The DNS takes about 100 s on 8 ranks.

## Expected outputs

- `sRe1cyl0.f00001` ... `sRe1cyl0.f00009` — leading SPOD modes (real part)
- `sIm1cyl0.f00001` ... `sIm1cyl0.f00009` — leading SPOD modes (imag part)
- `spod_stream_spectrum.dat` — SPOD eigenvalues (three blocks) in every frequency bin
- `mea1cyl0.f00001` — time-averaged (mean) flow field
- `logfile` — Slurm job log

## Result
Jean Zay, 8 ranks, 2026-10-09: the leading eigenvalue peaks in the bin St = 0.15625
(15.77); the bin St = 0.1875 holds 8.41. The shedding frequency (0.166) lies between the
bins, which are 0.03125 wide. These are the values of `../600_modal_pod`. `ref/reference.json`
holds them.

## How to run

```bash
mks 1cyl
sbatch run.local.slurm
```

The modal phase reads `1cyl0.f00001` through `1cyl0.f00100` from the working directory
(`modal_prefix='   '`, i.e. SESSION name, 100 snapshots at dt=0.5).
