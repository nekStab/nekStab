# Cylinder SFD Dynamic Testcase

This folder is the dynamic-tolerance SFD baseline for the cylinder base-flow
suite.

It is intentionally the non-OIFS case:

- `dyn` means tolerance scheduling through `ifdyntol = .true.`
- `oifs` is not active here
- the separate acceleration testcase lives in `../sfd_dyn_oifs`

## Physics

- 2D flow around a circular cylinder with `viscosity = -50.0`
  (`Re = 50` in Nek5000's negative-viscosity Reynolds-number convention)
- steady base flow computed with selective frequency damping

This testcase is intentionally Re = 50 so it can be compared directly with
`../sfd`.  Do not use this folder for the Re = 100 SFD experiment; the same
SFD parameters at Re = 100 do not give the documented monotone convergence.

## Active Case Definition

- `startFrom = BFRe40_1cyl0.f00001`
- shared seed hash:
  `77de81b1fa85b639a318546e36a2ef6049501c78fed284addc4a60259f48806f`
- `endTime = 400.0`
- `userParam01 = 1.1` for SFD
- `userParam04 = 0.12`
- `userParam05 = 0.05`
- `userParam08 = 0.9` for the dynamic velocity-tolerance factor
- `userParam09 = 100` for the dynamic tolerance update stride
- `variableDt = yes`
- `timeStepper = bdf3`
- `targetCFL = 0.5`
- `lx1 = 8`, `lxd = 12`
- `extrapolation` is left at the non-OIFS baseline
- pressure and velocity tolerances are `1e-9`
- `ew_tol_cap = 0` in `1cyl.usr`
- `viscosity = -50.0`

Dynamic tolerance scheduling is enabled in `1cyl.usr` with:

```fortran
ifdyntol = .true.
```

The cap is an optional upper bound.  Leave it at `0` for this benchmark so the
velocity tolerance follows `0.9 * residual` down to the final `1e-9`
convergence gate.  The scheduler updates every 100 timesteps and leaves the
pressure tolerance fixed.

## Latest Result

Latest local verification:

- date: 2026-04-28
- ranks: 8
- convergence time: `t = 168.8264`
- initial residual: `1.199387e-5` at `t = 0.005575878`
- final residual: `9.998440e-10`
- total timesteps: 30278
- elapsed wall time: 199.891 seconds
- scheduler cap: disabled (`ew_tol_cap = 0`)
- scheduler: velocity-only, `0.9 * residual`, update stride 100
- used tolerance range: `1e-9` to `4.991572e-6`
- output: converged `BF_1cyl0.f00001`

Conclusion: this case reaches the same strict gate as the fixed-tolerance
`../sfd` baseline while reducing measured 8-rank wall time from 231.736 seconds
to 199.891 seconds.  The scheduler path is active and logged; it relaxes only
the velocity Helmholtz tolerance and keeps the pressure tolerance fixed.

## Saved Artifacts

- `BFRe40_1cyl0.f00001`
  Initial condition used to start the testcase.
- `BFRe40_1cyl.nek5000`
  Collection file for visualization of the initial condition.
- `BF_1cyl0.f00001`
  Converged dynamic-tolerance SFD base flow.
- `BF_1cyl.nek5000`
  Collection file for visualization of the converged base flow.
- `plot.png`
  Current figure generated from the saved case data.
- `residu.png`
  Residual and scheduler-history figure.
- `residu.dat`
  SFD residual history used by the plots.
- `dyn_tol.dat`
  Scheduler history used to confirm the tolerance schedule.

## Plotting

To regenerate the figure:

```bash
python3 p_residu.py
python3 plot.py
```

The script writes:

- `plot.png`
- `residu.png` when `p_residu.py` is used
- `plot.png` overlays the fixed-tolerance `../sfd/residu.dat` baseline

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

- This folder should remain the clean dynamic-tolerance baseline.
- Do not enable `OIFS` here; that is a separate testcase.
- Keep both the initial condition and the actual saved result so the two states
  can be inspected or diffed directly.
- The seed restart time is reset to `t = 0` in `nekStab_usrchk` so residual and
  scheduler plots are referenced to the start of the testcase, not to the time
  stored in the Re = 40 seed file.
