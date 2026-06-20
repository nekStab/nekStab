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
`dm*cav0.f*` mode fields. The limit-cycle fundamental appears at **St = 0.0927**
(T = 10.785, matching the Floquet period) with `|μ| ≈ 1`, σ ≈ 0, plus harmonics
at 2×, 3×, … near the unit circle.

**Weakly-supercritical caveat:** Re = 1950 is only just above Re_c ≈ 1916, so
the cycle saturates very slowly (see [`../000_dns/README.md`](../000_dns/README.md)).
A weak residual transient (σ ≈ +0.004, `|μ| ≈ 1.002`) survives in the snapshots,
and because DMD ranks modes by norm, that transient mode is reported **above**
the St = 0.0927 fundamental. This is a property of the near-critical flow, not a
solver error — the fundamental and its harmonics are present lower in the
spectrum (the centre-probe FFT confirms St = 0.0927 is the true dominant
frequency). A textbook-clean ranking would need a much longer saturating DNS
(~3000 t.u.).

## Run
```bash
mks cav
sbatch run.local.slurm
```
