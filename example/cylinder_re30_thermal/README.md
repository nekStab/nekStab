# Thermally Coupled Cylinder at Re=30

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This case combines the circular-cylinder wake with scalar/thermal coupling and
Newton base-flow stages. The named thesis material discusses cylinder-wake
instability and thermosyphon scalar coupling separately, but it does not contain
a dedicated heated-cylinder case section. Thesis source for the cylinder
setting: chapter 2, `sec:ex:1cyl`; scalar-convection context:
`sec:ex:thermo`.

## Stages

- `000_dns_seed/` provides the DNS seed.
- `210_baseflow_newton/` computes the coupled fixed point.

## Temperature BCs from a velocity-only mesh (by design)

This case deliberately demonstrates a common real-world situation: the `.re2`
mesh carries **only velocity boundary conditions** (most meshes do). The
temperature BCs are not in the mesh — they are derived from the velocity BCs in
`usrdat2` (`cbc(:,:,2)` copied from `cbc(:,:,1)`, then walls and inflow set to
`'t '`).

For this to work the `.par` sets, under `[MESH]`, `numberOfBCFields = 1` so the
`.re2` reader reads just the one (velocity) BC field. Without it, enabling
temperature (`[TEMPERATURE]` → `ifheat`) makes the reader look for a second
(temperature) BC block the mesh does not contain, which aborts with
`Error reading .re2 boundary data ierr=8`. This is expected behaviour, not a
mesh defect — specifying temperature BCs in `usrdat2` is the point of the case.

