# Cylinder at Re=180

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This case targets the periodic cylinder wake below the 3D secondary threshold
and feeds Newton periodic-orbit, direct Floquet, adjoint Floquet, mode
animation, and OTD stages. The thesis uses the circular-cylinder wake to
connect periodic vortex shedding, Floquet multipliers, and pitchfork-type
secondary instability. Thesis source: chapter 2, `sec:ex:1cyl` and its
pitchfork-bifurcation subsection.

## Stages

- `000_dns_seed*/` provide periodic-wake seeds.
- `210_baseflow_newton/` computes the periodic orbit.
- `311_stability_direct_floquet/` and `321_stability_adjoint_floquet/` analyze
  perturbations over one shedding period.
- `410_postproc_animate_modes/` reconstructs mode dynamics.
- `500_otd/` tracks finite-time directions on the evolving wake.

