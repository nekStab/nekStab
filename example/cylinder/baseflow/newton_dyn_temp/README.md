# Cylinder Thermal Continuation Runbook

## Physics
2D flow around a circular cylinder with an active temperature field. The cylinder surface is heated (`T = 1` for `-1 < x < 1`), and buoyancy couples temperature to the y-momentum equation through the Boussinesq term in `userf`. `userParam06` is used as the Richardson/Rayleigh-like continuation parameter for the thermal coupling.

In this case:

- `Re` is the Reynolds number. It controls the ratio of inertia to viscosity. Higher `Re` makes the wake less damped and harder to converge with Newton.
- `Ri` is the thermal stratification / buoyancy-coupling coefficient. It multiplies the temperature field in the y-momentum forcing.
- `Ri = 0` recovers the isothermal cylinder case.
- `Ri > 0` means the heated fluid near the cylinder produces a positive upward buoyancy force. In an air-like interpretation, the cylinder heats the surrounding fluid, that warmer fluid rises, and the wake/plume is displaced upward.
- Larger `Ri` strengthens that upward thermal deflection, but it also makes the fixed-point solve more difficult.

This case is now used as a continuation workspace, not only as a one-shot Newton example. The stable path is:

1. Wash transients with DNS at low Reynolds number.
2. Use the washed thermal state as the Newton initial condition at the same Reynolds number.
3. Continue that converged thermal base flow to higher Reynolds number.
4. Only then run direct stability on the continued base flow.

Important interpretation:

- `Re = 20` is used here because the heated wake is still stable and therefore
  provides a practical Newton seed.
- `Re = 20` is not the final thermal-stability target if the goal is to study
  an unstable or nearly unstable heated wake.
- The direct modes of interest should be computed on the continued base flow
  (`Re = 25`, `Re = 30`, ...), not on the stable `Re = 20` seed itself.

## Critical Findings

- The temperature boundary conditions must be written to `cbc(:,:,2)` when `IFHEAT` is enabled. `ldimt` does not shift that slot.
- For restart files that must feed Newton, use `writeDoublePrecision = yes`. The old `#std 4` outputs caused restart trouble.
- Newton should use `variableDt = no`. The working setup recomputes a fixed `dt` from `endTime`.
- The Newton/Krylov norm already includes temperature. The runtime controls
  `thermal_norm_mode`, `thermal_buoyancy_coeff`, and `thermal_norm_weight`
  now support an automatic energy-based choice for the temperature weight.
  Pressure remains excluded.
- On this small 1996-element mesh, too many MPI ranks can corrupt the `elmap` in Newton checkpoint files (`nwt*` and `res*`). The symptom is the Pymech warning:
  `The 'elmap' appears to be corrupted as it contains an unexpected zero value.`
- This is not just a reader nuisance: ParaView can stop opening the files, and time is lost debugging fake "header" problems that are really bad parallel output for small local runs.
- Operational rule for local continuation tests on this mesh: keep the MPI count low. `8` ranks is a known-good choice here. Avoid `20` ranks for Newton checkpoint-producing runs unless the file layout has been revalidated.
- Opening files in ParaView can fail if a `.nek5000` collection file has stale `numtimesteps`. The raw `*.f0000*` files may still be fine.
- Starting Newton directly from a fully developed unstable thermal wake is much harder than starting from a washed low-Re thermal state.

## Step-by-Step Workflow

### Step 1. Thermal DNS Washout at Low `Re`

Use DNS to remove transients before Newton:

- `startFrom = 0`
- `userParam01 = 0`
- `userParam06 = 0.2` or `0.3`
- choose a Reynolds number low enough that the heated wake is still stable
- for the current continuation branch, start at `Re = 20`
- `viscosity = -20.0`
- `conductivity = -20.0`
- `variableDt = yes`
- `writeDoublePrecision = yes`
- run to `endTime = 100.0`

This produces a restartable thermal state:

- `1cyl0.f00001` at `t = 100`
- copy it to `BFt_1cyl0.f00001`
- run `visnek BFt_1cyl` so ParaView/VisIt can open the restart through
  `BFt_1cyl.nek5000`

### Step 2. Thermal Newton at the Same Low `Re`

Use the washed thermal DNS state as the Newton seed:

- `startFrom = BFt_1cyl0.f00001`
- `userParam01 = 2`
- `userParam06 = 0.2` or `0.3`
- keep the same low Reynolds number used for the washed DNS state
- `variableDt = no`
- use an Arnoldi sampling horizon tied to the wake frequency.
  Since `f ≈ 0.12` is only approximate, we round
  `1 / 0.12 / 8 ≈ 1.0417` to `endTime = 1.0`.
- large Krylov space: `userParam07 = 500`

#### Verified Result

At `Re = 20`, `Ri = 0.2`, the seeded thermal Newton solve converged:

- `NEWTON finished successfully after 6 iterations`
- output: `BF_1cyl0.f00001`

Preserve that converged base flow before continuing:

- copy `BF_1cyl0.f00001 -> BFre20t_1cyl0.f00001`
- run `visnek BF_1cyl`
- run `visnek BFre20t_1cyl`

#### Required Check Plot

After every Newton run, generate the convergence plot before judging the case:

```bash
uv run --with pymech --with matplotlib python plot.py
```

This figure is not optional in the workflow. It is the quickest way to inspect:

- `residu_newton.dat`
- `residu_arnoldi.dat`
- the converged base-flow velocity field
- the converged base-flow temperature field

### Step 3. Continue Toward `Re = 100`

Use the converged low-`Re` thermal Newton base flow as the restart and increase `Re` in small steps:

- `startFrom = BFre20t_1cyl0.f00001`
- `userParam01 = 2`
- keep `userParam06` fixed
- increase `Re` gradually, for example `20 -> 30 -> 40 -> 60 -> 80 -> 100`
- keep `variableDt = no`
- use `endTime = 1.0` as the practical Arnoldi sampling interval

#### Current Status

- `Re = 20`, `Ri = 0.2` Newton is now the first verified thermal base-flow continuation point.
- The next step is `Re = 30`, `Ri = 0.2` from `BFre20t_1cyl0.f00001`.
- That continuation step now uses the automatic clipped thermal weighting
  (`thermal_norm_mode = 2`) instead of a fixed manual value.
- Observed on 2026-04-15:
  the `Re = 30`, `Ri = 0.2` run starts cleanly and writes readable files on
  `8` ranks, but the first Newton linear solve stalls in the inner GMRES /
  Arnoldi loop. The residual decreases from `6.739824E-02` to about
  `5.664451E-02`, while the target for that linear solve is
  `1.684956E-02`.
- This means the automatic temperature weighting is not, by itself, enough to
  complete the `Re = 30` continuation step from the current `Re = 20` seed.
- Direct jumps to higher `Re` remain difficult and should still be avoided.

### Practical continuation rule

If `Re = 20` is known to be stable, use it only as the seed state. The next
useful workflow is:

1. `Re = 20` stable DNS washout
2. `Re = 20` Newton base flow
3. `Re = 25` Newton continuation
4. `Re = 30` Newton continuation
5. direct modes on the continued `Re = 30` base flow

This avoids confusing the stable seed case with the higher-`Re` case where the
mode structure is actually of interest.

## Files Used During Continuation

- `BFt_1cyl0.f00001`: washed thermal DNS seed
- `BFre20t_1cyl0.f00001`: converged `Re = 20` thermal Newton base flow
- `BF_1cyl0.f00001`: latest converged Newton base flow for the currently active Reynolds number

## Practical Notes

- Do not abort a Newton continuation just because the first residual goes up.
  In this thermal continuation path, the Newton / GMRES residual can increase
  before it eventually turns around and decreases.
- The correct stop rule is not "first increase = failure". Let the inner solve
  run long enough to show whether the residual is only non-monotone or has
  truly flattened at a plateau.
- If you change `thermal_norm_weight`, document the value in the case README or
  `.usr` because it changes the Newton/Krylov convergence metric.
- In this cylinder case, the buoyancy coefficient is `Ri`, because the forcing
  term is `Ri * T`. Automatic mode therefore uses
  `thermal_buoyancy_coeff = abs(uparam(6))`.
- With the current clipped automatic rule, the temperature field is still part
  of the Newton/Krylov norm, but the `Re = 30`, `Ri = 0.2` continuation run
  above suggests that a smaller continuation step or a different weighting rule
  is still needed there.
- The current implementation weights only the temperature field explicitly.
  Additional passive scalars still use unit weight for now, which keeps a clean
  extension path toward multi-scalar RANS cases later.
- If Pymech prints `The 'elmap' appears to be corrupted as it contains an unexpected zero value`, stop and suspect MPI count before changing headers or field-layout code.
- Always generate the `.nek5000` collection file for `BF*` outputs with
  `visnek`. The raw field file may be fine, but without the matching metadata
  file it is easy to think the restart is missing or unreadable in ParaView.
- In the tests here, the corruption showed up in `nwt*` and `res*` files at high rank count, while the main `BF_*` and DNS `1cyl*` files remained readable.
- If ParaView stops opening `nwt*` or `res*`, check the `.nek5000` collection file first before assuming the field files are corrupt.
- The `plot.py` script now tolerates empty `residu_newton.dat` files and will still render the available fields.
- The direct stability case under `example/cylinder/stability/direct` has already been patched to carry the temperature equation and plot the temperature component of the eigenvector.
