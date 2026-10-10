# Cylinder in forced transverse oscillation, Re = 100

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

The cylinder of diameter D moves across the stream with the displacement `A sin(2 pi f t)`, A = 1.08 D (`userParam05`) and
f = 0.1643 U/D (`userParam06`), at Re = 100. The mesh moves with the cylinder (`motion = user` in the `.par`): the `.usr`
sets the mesh velocity in `my_meshv`, and Nek5000 solves the flow in the arbitrary Lagrangian-Eulerian form. The forcing
frequency lies 1.2 % below the Strouhal number of the static cylinder on this mesh (0.1663, `000_dns`), and the wake locks to
the forcing: the flow is periodic with the period of the forcing, T = 1 / f = 6.0864 D/U.

Because the forcing sets the period, Newton-GMRES (`userParam01 = 2.2`) has no unknown period and no phase condition, and the
orbit has no phase mode: the Floquet multiplier 1 of an autonomous orbit is absent, and every multiplier lies inside the unit
circle when the lock-in is stable. The Floquet stage prints a warning that no multiplier lies near 1; the message says to ignore
it for a forced orbit.

The thesis chapters on periodic orbits and bifurcations of periodic orbits (chapter 2, `Periodic orbits`, `Bifurcations of
periodic orbits`, `sec:ex:1cyl`) give the concepts. They have no section on the forced oscillating cylinder; this case is
specific to the repository.

## Stages

| stage | mode | what it does |
|---|---|---|
| `000_dns` | 0 | the static cylinder (`userParam05 = 0`), with sponge and filter: the Strouhal number 0.1663 |
| `010_dns_forced` | 0 | DNS from rest with the forcing, 30 periods, one field per period: the wake locks to the forcing |
| `210_baseflow_newton` | 2.2 | Newton-GMRES from the last field of `010`: the forced orbit, 11 iterations to a squared residual of 1e-9 |
| `311_stability_direct_floquet` | 3.11 | direct Floquet analysis of that orbit: one converged multiplier, mu = 0.92988, real and inside the unit circle, so the lock-in is stable |

The stages run in this order, and each reads the field of the one before from its sibling folder (the field files are not
tracked). The orbit stages use the same `SIZE` (polynomial order 5, `lx1 = 6`, 2128 elements) and no filter, so that the DNS,
the Newton iteration and the Floquet analysis solve the same discrete equations.

## Moving mesh in the linearized solver

Nek5000 adds the term `B div(w u)` to the base flow when the mesh moves with the velocity `w`, and its perturbation solver
has no matching term. nekStab adds `B div(w u')` to the forcing of every perturbation (`ale_mesh_velocity_term` in
`src/forcing.f90`; DOC.md, section Moving Mesh). Without it the linearized map is not the derivative of the nonlinear
map: Newton on this orbit stopped at a squared residual of 1.9e-5, and the leading Floquet multiplier of a 30-vector check
was a pair of modulus 0.78 where finite differences of the nonlinear map gave 0.9319. With it Newton reaches 9e-10 in 11
iterations and the linearized and the finite-difference multipliers agree to 4.4e-4 over the eight leading ones
(`311_stability_direct_floquet/README.md`). The adjoint, transient-growth and OTD modes have no mesh-velocity term and stop
at start-up on a moving mesh.
