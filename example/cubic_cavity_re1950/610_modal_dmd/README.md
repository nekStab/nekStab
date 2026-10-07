# Cubic Cavity Re=1950 — Dynamic Mode Decomposition (DMD)

DMD of the saturated 3D limit cycle: extracts frequency-tagged modes and their
eigenvalues `μ = e^{(σ+iω)Δt}`.

## Snapshots
100 fields spanning the limit cycle, symlinked from the main DNS
(`cav0.f00001..100 -> ../000_dns/cav0.f00080..179`). Stored double precision
(single-precision 3D snapshots hang Nek MPI-IO on 16 ranks).

## nekStab Mode
`userParam01 = 6` — modal decomposition (reads the snapshot set, no time
stepping).

## Expected Output
`dmd_spectrum.dat` (`|μ|`, σ, ω, St, mode norms), `dmd_svd.dat`, and
`dm*cav0.f*` mode fields. The snapshots lie on the saturated limit cycle, so
every mode has `|μ| < 1` (σ between −4e-5 and −4.5e-3). The fundamental is at
**St = 0.09271** (T = 10.786; Floquet period 10.789) with `|μ| = 0.99981`. Its
harmonics at 2×, 3×, … are present with `|μ|` between 0.9998 and 0.9978. The
SVD holds 94.3 % of the energy in the first pair and 99.6 % in four modes.

DMD ranks the modes by norm. The 3rd harmonic (St = 0.278) has the largest
norm and ranks first; the fundamental ranks fifth. The centre-probe FFT gives
the same fundamental (St = 0.09269).

## Run
```bash
mks cav
sbatch run.local.slurm
```
