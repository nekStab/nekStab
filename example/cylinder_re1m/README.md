# Cylinder at Re=1e6

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This high-Reynolds-number cylinder directory is a finite-difference stability
scaffold for an open wake at a much higher operating point than the validated
low-Re cylinder examples. The named thesis chapters cover the circular-cylinder
wake as a canonical instability example but do not contain a dedicated Re=1e6
case section. Thesis source for concepts: chapter 2, `sec:ex:1cyl`; method
source: chapter 2, `Large-scale eigensolvers`.

## Stages

- `000_dns/` stores the nonlinear case scaffold.
- `310_stability_direct/findiff/` is the finite-difference direct-stability
  path for this operating point.

