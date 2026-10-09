# cylinder_re100 / 130_baseflow_dmt

Steady base flow of the Re = 100 cylinder wake by Dynamic Mode Tracking (DMT, Queguineur et al. 2019).
Compare with `../110_baseflow_sfd` (selective frequency damping) and `../120_baseflow_boostconv`.

## Method

DMT filters the velocity u with the equation `ū'' + β ū' + ω_c² ū = u'' + ω_c² u`. The filtered field ū is u without its content near the angular
frequency ω_c. The solver adds the force `f = −χ (u − ū)` to the momentum equation. The difference u − ū is the band-pass of u around ω_c, so the force damps an
oscillation at ω_c and does nothing to a steady flow. The band-pass has zero gain only at ω = 0, so a steady flow is the only state without force.
DMT therefore finds steady flows. It does not hold a periodic orbit: the orbit has content at its own frequencies, and the force acts on it.

## Case

- 2D flow around a circular cylinder, `viscosity = -100` (Re = 100), 2128 elements, `lx1 = 6`, `lxd = 9`.
- `startFrom = BF_seed_1cyl0.f00001`, the seed that the cylinder base-flow stages share. `1cyl.usr` sets the restart time to 0.
- `userParam01 = 1.3` (DMT). `userParam11` is ω_c = 0.7337, the frequency of the leading mode of the steady wake (stage `310_stability_direct`: σ = 0.1248,
  ω = 0.7337). `userParam12` is β = 1.0. `userParam13` is χ = 0.5, four times σ. `userParam14` = 5: the force starts at t = 5.
- The time step is constant (`variableDt = no`, `dt` from `targetCFL = 0.5`, about 0.01), because the filter uses the previous three steps. A change of the step
  stops the run with a message. Nek shortens the last step to land on `endTime`; DMT skips that step.
- The residual in `residu_dmt.dat` (column 2) is the root mean square over the domain of u − ū, divided by the time step. Column 3 is its change per step.
  The run stops when the residual falls below `userParam15 = 1e-7` and writes `BF_1cyl0.f00001` in double precision. `endTime = 800` is a limit that the run does
  not reach.

## Result (Jean Zay, 8 ranks, `ifort`)

- The residual meets the gate 1e-7 at t = 766.28 after 76228 time steps. It started at 2.1e-2. Wall time = 580 s.
- The converged field against the Newton solution of `../210_baseflow_newton/fp` (a different method for the same steady flow):

| statistic of the field | this run | Newton | relative difference |
|---|---|---|---|
| vx.min | -0.2001695 | -0.2001774 | 4.0e-05 |
| vx.max | 1.156066 | 1.156073 | 6.1e-06 |
| vx.norm | 264.3353 | 264.3375 | 8.3e-06 |
| vy.max | 0.6032385 | 0.6032416 | 5.1e-06 |

These are the differences of the SFD stages (4e-5, 6e-6, 8e-6, 5e-6).

## Reference

`ref/reference.json` checks six quantities. The last residual of `residu_dmt.dat` must lie below the gate 1e-7, so a run that stops at `endTime` fails. The time of the
last row must be within 5 % of 766.28. The four statistics of `BF_1cyl0.f00001` must be within 1e-4 of the Newton values above.
`python3 ../../../scripts/check_against_ref.py .` runs the check. The figures of the run are `plot_ic.png` (the seed), `plot_bf.png` (the converged field) and
`plot_residual.png`. `ref/` holds copies of them as `ic.png`, `bf.png` and `residual.png`.

## Run

```bash
mks 1cyl
sbatch run.local.slurm      # 8 ranks
python3 plot.py
python3 ../../../scripts/check_against_ref.py .
```

## Saved files

`BF_seed_1cyl0.f00001` is the initial condition. `BF_1cyl0.f00001`, `residu_dmt.dat` and the figures are regenerated.

## Literature

Queguineur et al. (2019), Phys. Fluids 31, 034101.
