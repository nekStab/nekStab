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
- `210_baseflow_newton/dyn_temp/` computes the coupled fixed point.

