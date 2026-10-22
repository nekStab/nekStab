# Cubic Cavity Re=1950 — Dynamic Mode Decomposition (DMD)

DMD of the saturated 3D limit cycle: extracts frequency-tagged modes and their
eigenvalues `μ = e^{(σ+iω)Δt}`.

## Snapshots
100 fields spanning the limit cycle, symlinked from the main DNS
(`cav0.f00001..100 -> ../000_dns/cav0.f00080..179`, t = 1740 to 1789.5, one every 0.5 time units, about 4.6 periods).
Stored double precision (single-precision 3D snapshots hang Nek MPI-IO on 16 ranks).

## nekStab Mode
`userParam01 = 6` — modal decomposition (reads the snapshot set, no time
stepping).

## Expected Output
`dmd_spectrum.dat` (`|μ|`, σ, ω, St, mode norms), `dmd_svd.dat`, and
`dm*cav0.f*` mode fields. The snapshots lie on the saturated limit cycle, so
every mode has `|μ| < 1` (σ between −4e-5 and −4.5e-3). The fundamental is at
**St = 0.09271** (T = 10.786; Floquet period 10.789) with `|μ| = 0.99981`. Its
harmonics n = 2 to 6 are present with `|μ|` between 0.99977 and 0.99998. The
SVD holds 94.3 % of the energy in the first two singular values and 99.6 % in four.

The list of modes is sorted by `|μ|`, that is by damping, not by frequency. The 3rd harmonic
(St = 0.278, `|μ|` = 0.999978) comes first; the fundamental is the fifth pair (modes 9 and 10). The centre-probe FFT gives
the same fundamental (St = 0.09269).

## Result
The run takes 11 s on 16 ranks of Jean Zay. `plot.py` takes, for each n = 1 to 6, the mode with the positive
frequency nearest to n 2π/T (T = 10.78909) and writes `dmd_harmonics.dat` and `plot_dmd.png`:

| n | mode in the list | ω / (2π/T) | σ (U/L) | norm |
|---|---|---|---|---|
| 1 | 9 | 1.00030 | −3.84e-4 | 0.779 |
| 2 | 7 | 1.99961 | −2.58e-4 | 0.716 |
| 3 | 1 | 3.00023 | −4.4e-5 | 0.716 |
| 4 | 13 | 4.00010 | −4.60e-4 | 0.705 |
| 5 | 3 | 4.99970 | −1.35e-4 | 0.708 |
| 6 | 5 | 6.00042 | −1.41e-4 | 0.705 |

Every one of the six modes lies on a harmonic of the cycle to 4.2e-4 or better. A cycle has σ = 0; the
fit over 4.6 periods leaves σ between −4.6e-4 and −4.4e-5. Modes 11 to 20 are the harmonics n = 7 and n = 8 (modes 11, 12
and 15, 16) and two pairs at ω / (2π/T) = 9.02 and 10.15 (modes 17 to 20). The last two pairs are not on a harmonic and
carry the strongest damping, σ = −3.6e-3 and −4.5e-3.

`ref/dmd.png`: (a) the norm of each mode against its frequency, with the harmonics marked, (b) the growth rate of each mode,
(c) to (e) the modulus of the velocity of the modes n = 1, 2, 3 on the plane z = 0.

## Reference
`ref/reference.json` checks twelve quantities: the largest `|μ|` (below 1 within 5e-5), the frequency of the six harmonic modes
(2e-3 of n 2π/T), the growth rate of the fundamental (5e-4), and the energy of each of the first four singular values
(0.1, 0.1, 0.05 and 0.05 percentage points). `python3 ../../../scripts/check_against_ref.py .` runs the check.

## Run
```bash
mks cav
sbatch run.local.slurm
uv run --with numpy --with matplotlib --with pymech python plot.py
```
