# Cylinder Re=100 — DMD modal analysis (uparam01=6.2)

## What this case shows

Dynamic Mode Decomposition extracts spatially coherent modes and their
associated frequencies and growth rates from a snapshot ensemble. Each DMD
mode has a single complex frequency; the leading mode at Re=100 captures the
vortex shedding at St ~= 0.164-0.167. Using `dmd_rank=20` avoids auto-rank
truncation that would miss the subdominant oscillatory modes.

## Bootstrap IC

`run.local.slurm` first runs a DNS (`dns.par`, `userParam01 = 0`) from the restart file
`rst_1cyl0.f00001` (a limit-cycle state at t = 150) to t = 200. It writes the 100 snapshots
`1cyl0.f00001` ... `1cyl0.f00100` (one every 0.5 time units), then restores `1cyl.par` and runs the
modal analysis on them. The restart file and `dns.par` are kept in `../600_modal_pod`; the script
copies them when they are missing. The DNS takes about 100 s on 8 ranks.

## Expected outputs

- `dm11cyl0.f00001`, `dm21cyl0.f00001` — leading DMD mode fields
- `dmd_spectrum.dat` — DMD eigenvalues (frequency, growth rate, amplitude)
- `dmd_svd.dat` — singular values of the snapshot matrix
- `mea1cyl0.f00001` — time-averaged (mean) flow field
- `logfile` — Slurm job log

## Result
Jean Zay, 8 ranks, 2026-10-09: the shedding pair has St = 0.16566 (second pair in the
norm ordering, after a pair near the Nyquist frequency), with harmonics at 0.33158 and
0.66301; the same values as in `../600_modal_pod`. `ref/reference.json` holds them.

## How to run

```bash
mks 1cyl
sbatch run.local.slurm
```

The modal phase reads `1cyl0.f00001` through `1cyl0.f00100` from the working directory
(`modal_prefix='   '`, i.e. SESSION name, 100 snapshots at dt=0.5).
