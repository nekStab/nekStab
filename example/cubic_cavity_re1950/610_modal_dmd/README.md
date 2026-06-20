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
`dm*cav0.f*` mode fields. For a limit cycle the leading DMD eigenvalues sit
**on the unit circle** (`|μ| ≈ 1.00`, σ ≈ 0) in conjugate pairs — the periodic
flow has no net growth/decay. The mode pairs recover the cycle's harmonics.

## Run
```bash
mks cav
sbatch run.local.slurm
```
