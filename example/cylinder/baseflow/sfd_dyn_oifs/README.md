# Cylinder SFD Dynamic OIFS Testcase

This folder tests selective frequency damping with dynamic solver tolerances
and OIFS time integration for the cylinder base-flow suite.

It is intentionally the accelerated comparison case:

- `dyn` means tolerance scheduling through `ifdyntol = .true.`
- `oifs` means the OIFS time-stepper acceleration is active
- OIFS is case-specific and must be validated against the non-OIFS baselines

## Physics

- 2D flow around a circular cylinder with `viscosity = -50.0`
  (`Re = 50` in Nek5000's negative-viscosity Reynolds-number convention)
- steady base flow computed with selective frequency damping

This testcase is intentionally Re = 50 so it can be compared directly with
`../sfd` and `../sfd_dyn`.  OIFS is an acceleration path, not a suite-wide
default.

## Active Case Definition

- `startFrom = BFRe40_1cyl0.f00001`
- shared seed hash:
  `77de81b1fa85b639a318546e36a2ef6049501c78fed284addc4a60259f48806f`
- `endTime = 400.0`
- `userParam01 = 1.1` for SFD
- `userParam04 = 0.12`
- `userParam05 = 0.05`
- `userParam08 = 0.95` for the dynamic velocity-tolerance factor
- `userParam09 = 100` for the dynamic tolerance update stride
- `variableDt = yes`
- `timeStepper = bdf3`
- `targetCFL = 5`
- `extrapolation = OIFS`
- `lx1 = 8`, `lxd = 12`
- pressure and velocity tolerances are `1e-9`
- `ew_tol_cap = 0` in `1cyl.usr`
- `viscosity = -50.0`

Dynamic tolerance scheduling is enabled in `1cyl.usr` with:

```fortran
ifdyntol = .true.
```

The cap is an upper bound.  Leave it at `0` for this benchmark so the velocity
tolerance follows `0.95 * residual` down to the final `1e-9` convergence gate.
The scheduler updates every 100 timesteps and leaves the pressure tolerance
fixed.

## Latest Result

Latest local verification:

- date: 2026-04-28
- ranks: 8
- convergence time: `t = 205.6384`
- initial residual: `9.863057e-5` at `t = 0.05575878`
- final residual: `9.982859e-10`
- total timesteps: 3688
- elapsed wall time: 110.790 seconds
- scheduler cap: disabled (`ew_tol_cap = 0`)
- scheduler: velocity-only, `0.95 * residual`, update stride 100
- used tolerance range: `1.294030e-9` to `1.753468e-5`
- output: converged `BF_1cyl0.f00001`

Comparison against the non-OIFS baselines from the same seed:

- `../sfd`: 30268 timesteps, final residual `9.996892e-10`,
  elapsed wall time 231.736 seconds
- `../sfd_dyn`: 30278 timesteps, final residual `9.998440e-10`,
  elapsed wall time 199.891 seconds
- `../sfd_dyn_oifs`: 3688 timesteps, final residual `9.982859e-10`,
  elapsed wall time 110.790 seconds

Conclusion: OIFS reaches the same strict `1e-9` convergence gate with about
8.2 times fewer timesteps and about 1.8 times less wall time than `../sfd_dyn`
for this Re = 50 testcase.

The step-count speedup and the wall-clock speedup are intentionally reported
separately.  OIFS advances much farther per timestep, but each OIFS timestep is
more expensive than the standard BDF timestep.  The plotting scripts therefore
include a measured-cost panel using the verified 8-rank wall times above:
`../sfd_dyn_oifs` costs about `0.55` times the wall time of `../sfd_dyn`, not
`1 / 8.2` times.  If the cases are rerun, update the elapsed-time constants in
`plot.py` and `p_residu.py` together with this section.

## Saved Artifacts

- `BFRe40_1cyl0.f00001`
  Initial condition used to start the testcase.
- `BFRe40_1cyl.nek5000`
  Collection file for visualization of the initial condition.
- `BF_1cyl0.f00001`
  Converged dynamic-tolerance OIFS SFD base flow.
- `BF_1cyl.nek5000`
  Collection file for visualization of the converged base flow.
- `plot.png`
  Current figure generated from the saved case data.
- `residu.png`
  Residual comparison, measured-cost comparison, and scheduler-history figure.
- `residu.dat`
  SFD residual history used by the plots.
- `dyn_tol.dat`
  Scheduler history used to confirm the active tolerance cap.

## Running And Plotting

Run the testcase from this folder with 8 ranks:

```bash
mks 1cyl
mpiexec -np 8 ./nek5000
visnek BFRe40_1cyl
visnek BF_1cyl
python3 p_residu.py
python3 plot.py
```

The scripts write:

- `residu.png`, which overlays `../sfd`, `../sfd_dyn`, and this OIFS case and
  shows both physical-time convergence, measured 8-rank wall cost, and the
  requested/used solver tolerances from `dyn_tol.dat`
- `plot.png`, which shows the physical-time residual comparison, measured
  8-rank wall-cost comparison, and the saved final field

The dynamic scheduler also writes:

- `dyn_tol.dat`

Its columns are:

- time
- SFD residual
- current Nek5000 solver tolerance before the update
- requested scheduler tolerance `userParam08 * residual`
- tolerance actually used after applying any cap
- cap value

## Notes

- Keep this folder separate from `../sfd_dyn`; OIFS must not be enabled in the
  strict dynamic-tolerance baseline.
- `SESSION.NAME` must point to this folder.  If it points to `../sfd_dyn`, the
  run will silently read and overwrite the wrong case artifacts.
- Keep both the initial condition and the actual saved result so the two states
  can be inspected or diffed directly.
- The seed restart time is reset to `t = 0` in `nekStab_usrchk` so residual and
  scheduler plots are referenced to the start of the testcase, not to the time
  stored in the Re = 40 seed file.
