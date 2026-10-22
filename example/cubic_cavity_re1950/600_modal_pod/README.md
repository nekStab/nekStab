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
concentrates in the **leading real+imag pair**: mode 1 = 57.6 %, mode 2 = 36.7 %
(94.3 % cumulative by two modes) — the oscillation captured as one POD pair.

## Result
The run takes 8 s on 16 ranks of Jean Zay. `plot.py` reads `pod_spectrum.dat`, `pod_coefficients.dat` and the two leading modes,
and writes `pod_summary.dat` and `plot_pod.png`. The energy of the first ten modes comes in pairs:

| modes | energy (%) | cumulative (%) | turning rate of the pair, ω / (2π/T) |
|---|---|---|---|
| 1, 2 | 57.558, 36.708 | 94.266 | 1.0005 |
| 3, 4 | 2.779, 2.559 | 99.605 | 1.9998 |
| 5, 6 | 0.183, 0.172 | 99.959 | 2.9989 |
| 7, 8 | 0.0184, 0.0181 | 99.996 | 3.9981 |
| 9, 10 | 0.0018, 0.0018 | 99.9996 | 4.9975 |

The turning rate is the rate of the unwrapped phase atan2(a₂, a₁) of the coefficients of the pair, fitted over the 100 snapshots
(4.6 periods, T = 10.78909). The pairs turn at the harmonics 1 to 5 of the cycle, so the POD of this flow separates
the harmonics in time. The two modes of the first pair do not have the same energy because the orbit of (a₁, a₂) is an ellipse.

`ref/pod.png`: (a) the energy of the first twelve modes, (b) the orbit of the coefficients (a₁, a₂), (c) and (d) the speed of modes 1 and 2
on the plane z = 0.

## Reference
`ref/reference.json` checks twelve quantities: the energy of modes 1 to 4 (0.05 and 0.02 percentage points), the energy of the first pair and of the
first four modes, the turning rate of the five pairs (2e-3 for the first, 3e-3 for the others), and the largest eigenvalue (0.5 %).
`python3 ../../../scripts/check_against_ref.py .` runs the check.

## Run
```bash
mks cav
sbatch run.local.slurm
uv run --with numpy --with matplotlib --with pymech python plot.py
```
