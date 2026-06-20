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

- `startFrom` is commented out: the archived `BFre25t_1cyl0.f00001` seed
  belonged to the legacy `nelg = 1996` mesh and is incompatible with the
  current `nelg = 2128` unified cylinder mesh (see `NOTE-mesh-replaced-*.md`).
  The case now cold-starts from `useric` (`ux = 1`, `T = 0`).
- `userParam01 = 2` for Newton-GMRES
- `userParam06 = 0.2`
- `userParam07 = 200`
- `endTime = 1.1`
- `variableDt = yes`
- pressure, velocity, and temperature tolerances are all `1e-8`
- `viscosity = -30.0`
- `conductivity = -30.0`

The archived seed and previous solved outputs remain in the folder as
references:

- `BFre25t_1cyl0.f00001` and `BFre25t_1cyl.nek5000`
  Legacy seed thermal base flow at `Re = 25`, mesh `nelg = 1996`. Not loaded
  by the current `.par`.
- `BF_1cyl0.f00001` and `BF_1cyl.nek5000`
  Solved thermal base flow at `Re = 30`, produced before the mesh
  unification. Kept for plot regeneration.

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

This case is **also a worked example of adding a temperature field to a mesh
that does not carry temperature boundary conditions**. The shipped `1cyl.re2`
is a standard non-thermal cylinder mesh — `genmap` only writes the velocity
slice of `cbc`. Temperature BCs are patched in entirely from `usrdat2`:

```fortran
if (ifheat) then
  cbc(:,:,2) = cbc(:,:,1)            ! clone velocity BCs as starting point
  do iel = 1,nelt
    do ifc = 1, 2*ndim
      if (cbc(ifc,iel,1).eq.'W  ') cbc(ifc,iel,2) = 't  '   ! no-slip wall  -> Dirichlet temp
      if (cbc(ifc,iel,1).eq.'v  ') cbc(ifc,iel,2) = 't  '   ! velocity inflow -> Dirichlet temp
    enddo
  enddo
endif
```

This is intentional and demonstrates the recommended pattern when you need
to add a scalar field on top of a velocity-only mesh: edit `cbc(:,:,2)` in
`usrdat2`, then supply the actual values in `userbc`. Here the cylinder
surface receives `T = 1` on `-1 < x < 1` and `T = 0` elsewhere; buoyancy
enters y-momentum as `Ri * T` in `userf`.

Key points:

- The temperature BC slot is always `cbc(:,:,2)` when `IFHEAT` is enabled.
  Increasing `ldimt` adds scalar slots at indices `3, 4, ...` and does
  **not** move temperature off slot `2`.
- The internal BC codes `'E  '` and `'P  '` (interior face / periodic) are
  inherited from the velocity slice by the initial `cbc(:,:,2) = cbc(:,:,1)`
  copy, so the loop only needs to handle physical boundaries.
- Without the `usrdat2` patch, a thermal run on this mesh would fail with
  `ierr=8` ("Error reading .re2 boundary data") because Nek5000 looks for
  per-face temperature codes that are not in the file.

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
