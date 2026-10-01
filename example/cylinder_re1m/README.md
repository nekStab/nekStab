# Cylinder at Re=1e6

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This high-Reynolds-number cylinder directory carries the RANS stability
pipeline at a much higher operating point than the validated low-Re cylinder
examples. The named thesis chapters cover the circular-cylinder wake as a
canonical instability example but do not contain a dedicated Re=1e6 case
section. Thesis source for concepts: chapter 2, `sec:ex:1cyl`; method source:
chapter 2, `Large-scale eigensolvers`.

## Stages

- `000_dns/` is the nonlinear k-tau march on the 2D mesh (`lelg=1480`).
- `110_baseflow_sfd/` builds the SFD mean the stability stages restart from.
- `310_stability_direct/coupled/` is the coupled operator (full RANS Fréchet).
- `310_stability_direct/quasilaminar/` is the quasilaminar operator (base eddy viscosity prescribed).

The run mesh is each stage's `1cyl.re2` (1480 elements, body-fitted, closest node at r = 0.5). The generator `310_stability_direct/quasilaminar/mesh/1cyl2d.re2` (1464 elements) is Cartesian through the disk, with nodes at r = 0 and a 0.3 pitch. Same outer domain. Regenerating from that generator does not reproduce the run mesh.
