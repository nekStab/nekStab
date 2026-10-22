# Cylinder Re=100 — DMD modal analysis (uparam01=6.2)

## What this case shows

Dynamic Mode Decomposition extracts spatially coherent modes and their
associated frequencies and growth rates from a snapshot ensemble. Each DMD
mode has a single complex frequency. At Re=100 the snapshots lie on the limit cycle, so
the eigenvalues sit on the shedding frequency St = 0.1657 and on its multiples.
Using `dmd_rank=20` avoids auto-rank truncation that would miss the subdominant oscillatory modes.

## Bootstrap IC

`run.local.slurm` first runs a DNS (`dns.par`, `userParam01 = 0`) from the restart file
`rst_1cyl0.f00001` (a limit-cycle state at t = 150) to t = 200. It writes the 100 snapshots
`1cyl0.f00001` ... `1cyl0.f00100` (one every 0.5 time units), then restores `1cyl.par` and runs the
modal analysis on them. The restart file and `dns.par` are kept in `../600_modal_pod`; the script
copies them when they are missing. The DNS takes about 100 s on 8 ranks.

## Expected outputs

- `dm11cyl0.f00001` to `dm11cyl0.f00003` and `dm21cyl0.f00001` to `dm21cyl0.f00003` — real and imaginary
  part of the modes 1 to 3 of `dmd_spectrum.dat` (`modal_nsave = 3` in `1cyl.usr`). Mode 3 is the shedding mode
- `dmd_spectrum.dat` — DMD eigenvalues (frequency, growth rate, amplitude)
- `dmd_svd.dat` — singular values of the snapshot matrix
- `mea1cyl0.f00001` — time-averaged (mean) flow field
- `logfile` — Slurm job log

## Result
Jean Zay, 8 ranks, 2026-10-09: `dmd_spectrum.dat` lists the 20 eigenvalues of rank 20, ordered by |mu|. The first ten
are five conjugate pairs with |mu| from 0.99982 to 0.99917 and growth rates sigma between -3.7e-4 and -1.7e-3 per time
unit. Their Strouhal numbers are 0.49736, 0.16566, 0.33158, 0.66301 and 0.83965. The first four are 3.002, 1, 2.002 and
4.002 times the shedding frequency 0.16566; the values 0.16566, 0.33158 and 0.66301 are those of `../600_modal_pod`.
The other ten eigenvalues are damped more strongly (sigma from -1.9e-3 to -0.23). `ref/reference.json` holds the
frequencies of the first four pairs.

`plot.py` writes `plot_spectrum.png` (growth rate against Strouhal number of the eigenvalues, with the multiples of the
shedding frequency), `plot_svd.png` (singular values of the snapshot matrix, with the rank used) and `plot_mode.png` (vy of
the real part of the shedding mode). `ref/` holds copies as `spectrum.png`, `svd.png` and `mode.png`.

## How to run

```bash
mks 1cyl
sbatch run.local.slurm
```

The modal phase reads `1cyl0.f00001` through `1cyl0.f00100` from the working directory
(`modal_prefix='   '`, i.e. SESSION name, 100 snapshots at dt=0.5).
