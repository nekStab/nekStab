# Cylinder Thermal Newton Testcase

This folder is a clean Newton testcase for the heated-cylinder problem with an
active temperature field and dynamic inner tolerances.

The current saved jump is:

- seed base flow: `Re = 25`, `Ri = 0.2`
- target base flow: `Re = 30`, `Ri = 0.2`

This folder is only for the steady thermal Newton solve. It is not a direct
stability case and it is not a multi-step continuation ladder.

## Physics

- `Re` is the Reynolds number.
- `Ri` is the thermal buoyancy coupling used in `userf`.
- `Ri > 0` means the cylinder heats the fluid and the warmer fluid rises,
  which shifts the wake upward.
- `Ri = 0` recovers the isothermal cylinder case.

## Active Case Definition

The current `.par` file is configured as follows:

- `startFrom = BFre25t_1cyl0.f00001`
- `userParam01 = 2` for Newton-GMRES
- `userParam06 = 0.2`
- `userParam07 = 200`
- `endTime = 1.1`
- `variableDt = yes`
- pressure, velocity, and temperature tolerances are all `1e-8`
- `viscosity = -30.0`
- `conductivity = -30.0`

The current restart chain in this folder is:

- `BFre25t_1cyl0.f00001` and `BFre25t_1cyl.nek5000`
  Seed thermal base flow at `Re = 25` saved at `t = 1.0`.
- `BF_1cyl0.f00001` and `BF_1cyl.nek5000`
  Solved thermal base flow at `Re = 30` saved at `t = 1.1`.

## Files In This Folder

Case-defining files:

- `1cyl.par`
- `1cyl.usr`
- `1cyl.re2`
- `1cyl.ma2`
- `SIZE`

Saved solution files:

- `BFre25t_1cyl0.f00001`
- `BFre25t_1cyl.nek5000`
- `BF_1cyl0.f00001`
- `BF_1cyl.nek5000`

Plotting files:

- `plot.py`
- `plot_re25_to_re30_convergence.png`

Generated helper state:

- `.state`

The archived convergence figure is kept in the folder. The raw
`residu_newton.dat`, `residu_gmres.dat`, and `residu_arnoldi.dat` files are
not part of the current clean payload.

## What `1cyl.usr` Is Testing

This folder is named `newton_dyn_temp` because all three parts are active:

- `newton`: steady state is computed with Newton-GMRES
- `dyn`: `ifdyntol = .true.` is enabled
- `temp`: temperature is part of the coupled solve

The user file also contains the current thermal-norm setup:

- `thermal_norm_mode = 2`
- `thermal_buoyancy_coeff = abs(uparam(6))`
- `thermal_norm_min = 1.0`
- `thermal_norm_max = 10.0`

This keeps temperature inside the Newton/Krylov norm while preventing the
thermal contribution from becoming excessively stiff in this cylinder case.

## Boundary Conditions And Coupling

- The cylinder surface is heated with `T = 1` on `-1 < x < 1`.
- Buoyancy is added in the y-momentum equation through `Ri * T`.
- Temperature boundary conditions are written to `cbc(:,:,2)` when `IFHEAT` is
  enabled.
- Increasing `ldimt` does not move that temperature boundary-condition slot.

## Plotting

To regenerate the current field figure:

```bash
uv run --with pymech --with matplotlib python plot.py
```

The script writes:

- `plot.png`

The current figure style is intentional:

- panel `(a)` is the convergence panel when residual logs are present
- velocity uses a divergent colormap centered at `|u| = 1`, so the free stream
  is white
- temperature uses a white-to-warm colormap, so `T = 0` is white
- the two scalar colorbars are outside the axes and placed on top

If no residual logs are present, `plot.py` still plots the saved base-flow
fields from the latest available `BF*` file.

## Practical Notes

- Keep the MPI count low on this small mesh. Too many ranks can corrupt the
  `elmap` in Newton-related output files.
- If Pymech warns that `elmap` contains an unexpected zero value, suspect the
  MPI partition before changing file headers or restart logic.
- Always run `visnek` for the saved `BF*` outputs so ParaView or VisIt can
  open the matching `.nek5000` collection file.
- Do not treat a first residual increase as immediate failure. In this case
  the Newton residual can increase before it decreases.
