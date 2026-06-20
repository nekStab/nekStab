# Cubic Cavity Re=1950 — Proper Orthogonal Decomposition (POD)

POD of the saturated 3D limit cycle: extracts the energy-ranked coherent
structures of the post-Hopf flow.

## Snapshots
100 fields spanning the limit cycle, symlinked from the main DNS
(`cav0.f00001..100 -> ../000_dns/cav0.f00080..179`). Stored double precision
(single-precision 3D snapshots hang Nek MPI-IO on 16 ranks).

## nekStab Mode
`userParam01 = 6` — modal decomposition (reads the snapshot set, no time
stepping).

## Expected Output
`pod_spectrum.dat` (eigenvalue / energy / cumulative), `pod_coefficients.dat`,
and `pod*cav0.f*` mode fields. For a single-frequency limit cycle the energy
concentrates in the **leading real+imag pair**: mode 1 ≈ 58 %, mode 2 ≈ 36 %
(≈ 94 % cumulative by two modes) — the oscillation captured as one POD pair.

## Run
```bash
mks cav
sbatch run.local.slurm
```
