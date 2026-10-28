# cylinder_re100 / 110_baseflow_sfd / casacuberta_dyn

Steady base flow of the Re = 100 cylinder wake by selective frequency damping, Casacuberta formulation, dynamic solver tolerance.
Compare with ../casacuberta (fixed tolerance) and ../casacuberta_dyn_oifs (OIFS); the family README (`../README.md`) has the table of all six.

## Case

- 2D flow around a circular cylinder, `viscosity = -100` (Re = 100), 2128 elements, `lx1 = 6`, `lxd = 9`.
- `startFrom = BF_seed_1cyl0.f00001`, the seed that the cylinder base-flow stages share. `1cyl.usr` sets the restart time to 0.
- `userParam01 = 1.1` (SFD), `userParam04 = -0.166` (negative: Casacuberta), `userParam05 = 0.26`.
  `userParam04` is the frequency in cycles per time unit (omega = 2 pi f = 1.043) and `userParam05` the growth rate sigma
  of the formulas in `DOC.md`. The run prints `cutoff 0.4075, gain -0.6675`.
- The time stepper is BDF3 with `targetCFL = 0.5` (step about 0.01, variable).
- `ifdyntol = .true.` in `1cyl.usr`: the velocity solver tolerance follows `userParam08` = 0.9 times the residual, is updated every `userParam09` = 100 steps and never falls below the gate; the pressure tolerance stays fixed. `dyn_tol.dat` records the schedule (time, residual, tolerance before the update, requested tolerance, tolerance used, cap).
- The run stops when the residual falls below the larger `residualTol` of the `.par` (1e-10) and writes `BF_1cyl0.f00001` in double
  precision. `endTime = 1500` is a limit that the run does not reach.
- The residual in `residu.dat` (column 2) is the root mean square over the domain of the change of the velocity during one time step.

## Result (Jean Zay, 8 ranks, `ifort`)

- The residual meets the gate 1e-10 at t = 869.16 after 86462 time steps; it ends at 9.996e-11
  and started at 2.09e-04. Wall time 588 s.
- Over the last fifth of the run the residual falls as exp(-0.0177 t).
- The converged field against the Newton solution of `../../210_baseflow_newton` (a different method for the same steady flow):

| statistic of the field | this run | Newton | relative difference |
|---|---|---|---|
| vx.min | -0.2001695 | -0.2001774 | 4.0e-05 |
| vx.max | 1.156066 | 1.156073 | 6.1e-06 |
| vx.norm | 264.3354 | 264.3375 | 8.1e-06 |
| vy.max | 0.6032385 | 0.6032416 | 5.2e-06 |

The strict gate brings the field to 4e-5 of the Newton solution in every statistic; the tolerance of the reference (1e-4) covers that.

## Reference

`ref/reference.json` checks six quantities. The last residual of `residu.dat` must lie below the gate (tolerance equal to the gate, so a
run that stops at `endTime` fails). The time of the last row must be within 5 % of 869.16. The four statistics of `BF_1cyl0.f00001`
must be within 1e-4 of the Newton values above. `python3 ../../../../scripts/check_against_ref.py .` runs the check; the figures of the
run are `plot_ic.png` (the seed), `plot_bf.png` (the converged field), `plot_residual.png`, `plot_dynamic_tolerance.png`. `ref/compare.png` (made by `../plot_compare.py` from the runs of 2026-10-28) overlays the residual of the fixed, dyn and dyn OIFS runs against the time step and against the wall time.

## Run

```bash
mks 1cyl
sbatch run.local.slurm      # 8 ranks
python3 plot.py
python3 ../../../../scripts/check_against_ref.py .
```

Use 8 ranks. Higher rank counts on this small mesh can write field files with corrupted `elmap` metadata.

## Saved files

`BF_seed_1cyl0.f00001` is the initial condition. `BF_1cyl0.f00001` is a converged base flow kept as proof.
`residu.dat` and `dyn_tol.dat` are regenerated.

## Literature

Åkervik et al. (2006), Phys. Fluids 18, 068102; Casacuberta et al. (2018), J. Comput. Phys. 375.
