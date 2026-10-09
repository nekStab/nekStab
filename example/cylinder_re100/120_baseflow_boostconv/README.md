# cylinder_re100 / 120_baseflow_boostconv

Steady base flow of the cylinder wake at Re = 50 by BoostConv (`userParam01 = 1.2`, Citro et al. 2017). The folder belongs to the
`cylinder_re100` family, but this stage runs at Re = 50 (section "Why Re = 50"). Compare with `../110_baseflow_sfd` and `../130_baseflow_dmt`.

## Method

Plain time marching does not reach an unstable steady flow. BoostConv takes the change of the velocity over one time step as the residual.
Every `bst_skp` steps it stores the residual. From the last `bst_snp` stored residuals it forms differences and orthonormalises them (modified
Gram-Schmidt, weighted with the mass matrix). It finds the combination of these differences that cancels the new residual in the least-squares
sense, and corrects the velocity update by the same combination of the stored updates. This is an Anderson-type extrapolation. A steady flow
has zero residual, so it is the fixed point.

BoostConv acts on the velocity only. With a temperature or another scalar the code prints a warning and does not accelerate the scalar; use SFD
for such a case.

## Case

- 2D flow around a circular cylinder, `viscosity = -50` (Re = 50), 2128 elements, `lx1 = 6`, `lxd = 9`, explicit filter of weight 0.01, 8 ranks (`lpmin = 8`).
  The critical Reynolds number of the wake is 46.6, so the steady flow at Re = 50 is unstable.
- `startFrom = BF_seed_1cyl0.f00001`, the seed that the cylinder base-flow stages share: a snapshot of the shedding wake at Re = 100 (`../000_dns`).
- `userParam01 = 1.2`. The subspace size `bst_snp` and the sampling interval `bst_skp` are not in the `.par`. The defaults of `src/main.f90`
  (10 and 10) apply here. To change them, set both in `nekStab_usrchk` of `1cyl.usr`.
- The residual in `residu.dat` (column 2) is the root mean square over the domain of the change of the velocity during one time step, sampled
  every `bst_skp` steps. The run stops when it falls below the larger of the two `residualTol` of the `.par` (1e-6) and writes `BF_1cyl0.f00001` in double
  precision. `endTime = 500` is a limit.

## Result (Jean Zay, 8 ranks, `ifort`, 2026-10-09)

The residual meets the gate 1e-6 at t = 159.43, after 15 860 time steps (1586 samples). It started at 2.1e-4. Wall time 146 s (150 s and 153 s in two repeats; the
repeats gave the same `residu.dat` to every digit). The converged field against the Newton solution at Re = 50 (`../210_baseflow_newton/fp` with
`viscosity = -50`: five iterations, residual 9.6e-10):

| statistic of the field | this run | Newton | relative difference |
|---|---|---|---|
| vx.min | -0.1220497 | -0.1220482 | 1.2e-05 |
| vx.max | 1.164340 | 1.164336 | 3.6e-06 |
| vx.norm | 262.6001 | 262.5338 | 2.5e-04 |
| vy.max | 0.5636222 | 0.563611 | 2.0e-05 |

The norm sums the whole domain. The difference of 2.5e-4 is the part of the field that the gate 1e-6 leaves unconverged: with the gate 1e-8 the
norm differs by 8.0e-7 (table below).

## Parameters

Same case and seed, Re = 50, Jean Zay, 8 ranks. The gate is 1e-6 except in the last two rows. "Difference" is the relative difference of `vx.norm`
from the Newton value.

| `bst_snp` | `bst_skp` | gate | converged at t | time steps | wall time (s) | difference |
|---|---|---|---|---|---|---|
| 10 (default) | 10 (default) | 1e-6 | 159.4 | 15 860 | 146 | 2.5e-04 |
| 10 | 5 | 1e-6 | 152.4 | 15 160 | 156 | 6.8e-04 |
| 30 | 5 | 1e-6 | 168.3 | 16 745 | 221 | 6.8e-04 |
| 30 | 10 | 1e-6 | not reached by t = 500 | 49 730 | 535 and 634 | |
| 10 (default) | 10 (default) | 1e-8 | 232.4 | 24 420 | 310 | 8.0e-07 |
| 30 | 5 | 1e-8 | 292.0 | 29 050 | 649 | 3.4e-06 |

- Three settings reach 1e-6 between t = 152 and 168. A subspace of 30 is not better than one of 10: at `bst_skp = 5` it needs t = 168 against 152 and
  221 s against 156 s.
- With `bst_snp = 30` and `bst_skp = 10` the run does not reach the gate. The residual at t = 500 is 4.1e-6 in one run and 7.9e-6 in a repeat. Both
  runs have the same source and input, ran on different nodes, and have residuals that differ from time step 240 on. Repeats of the other settings
  (10 and 10, 10 and 5, 30 and 5) gave the same time of convergence to every digit.
- The gate sets the accuracy of the field. Going from 1e-6 to 1e-8 costs 46 % more time for the default and 73 % more for 30 and 5. The norm then
  agrees with Newton to 8.0e-7 and 3.4e-6.
- Wall times of single runs vary by tens of percent on Jean Zay (the default needed 9 ms per time step in the gate run and 13 ms in the run with the
  gate 1e-8). Compare the number of time steps, not the wall time.

## Why Re = 50

At Re = 100, the Reynolds number of the family, BoostConv did not reach the gate in this setup. With the default subspace (10 and 10) and with 30 and 5, the
residual stayed between 7e-5 and 9e-4 up to t = 500. Its median over the last 100 time units was 2.3e-4 for the default and 1.9e-4 for 30 and 5. The stage
therefore runs at Re = 50, just above the onset of shedding. For the steady wake at Re = 100 use `../110_baseflow_sfd`, `../130_baseflow_dmt` or
`../210_baseflow_newton`.

## Reference and figures

`ref/reference.json` checks six quantities. The last residual of `residu.dat` must lie below the gate 1e-6, so a run that stops at `endTime` fails.
The time of the last row must be within 25 % of 159.43. The four statistics of `BF_1cyl0.f00001` must be within 2e-4 (extrema) and 3e-3 (norm) of
the Newton values above. `python3 ../../../scripts/check_against_ref.py .` runs the check. The figures of the run are `plot_ic.png` (the seed),
`plot_bf.png` (the converged field) and `plot_residual.png` (the residual against time, with the gate). `ref/` holds copies of them as `ic.png`,
`bf.png` and `residual.png`.

## Run

```bash
mks 1cyl
sbatch run.local.slurm      # 8 ranks
python3 plot.py
python3 ../../../scripts/check_against_ref.py .
```

## Saved files

`BF_seed_1cyl0.f00001` is the initial condition. `BF_1cyl0.f00001`, `residu.dat` and the figures are regenerated.

## Literature

Citro, Luchini, Giannetti, Auteri (2017), J. Comput. Phys. 344, 234-246.
