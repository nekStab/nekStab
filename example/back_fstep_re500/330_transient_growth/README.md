# back_fstep_re500/330_transient_growth

Optimal transient growth of the backward-facing step at Re=500
(`userParam01 = 3.3`). The base flow is stable, but the shear layer behind
the step amplifies a well-chosen initial perturbation by more than four
orders of magnitude before the packet leaves the domain. For a horizon tau,
nekStab computes the largest energy gain

    G(tau) = max over q0 of  E(tau) / E(0)

as the leading eigenvalue of M^H M, where M = exp(tau L) is the linearized
flow map and M^H its adjoint. Each Krylov step is one forward and one
adjoint integration over tau.

## Setup

- Base flow: `BF_bfs0.f00001`, the Newton solution of
  `../210_baseflow_newton/` (tracked here, so this stage runs on its own).
- `endTime = 57.905` is the horizon tau. This value is the peak of the
  reference envelope.
- `userParam07 = 64`: Krylov basis size. `dt = 0.005` is fixed, so the
  forward and adjoint integrations cover the same tau.
- Sponge: 5 units at the inlet, 5 units at the outlet (`userParam08/09`).
  A 10-unit outlet sponge damps the packet before tau = 58.
- 8 MPI ranks (`lpmin = 8`). Wall time at tau = 57.905: 8 h 40 min in one
  run, about 20 h in another and 23 h in a third, all on 8 ranks of a 6-core
  node (the speed differs between runs). The third run reached a 23 h limit
  a few minutes after it wrote the results. `run.local.slurm` asks for 24 h.

Run: `sbatch run.local.slurm` (or `mpiexec -np 8 ./nek5000`).
G is the first value in `Spectre_Hp_conv.dat`. Plot: `python plot.py`.

## Result

One run per horizon, at the reference abscissae of Blackburn, Barkley &
Sherwin (2008), fig. 5 (`barkley2008_fig5.ref`). Values in
`growth_sweep.csv`:

| tau | G nekStab | G reference | difference |
|---:|---:|---:|---:|
| 19.916 | 2535 | 1997 | +27 % |
| 29.913 | 13127 | 12273 | +7.0 % |
| 39.838 | 30605 | 32250 | -5.1 % |
| 57.905 | 61690 | 63152 | -2.3 % |

The peak gain agrees to 2.3 %. The reference values come from a curve read
off a log-scale figure. At tau = 20 the curve rises by a factor 1.6 per
2 time units, so the +27 % difference equals a shift of about 0.8 in tau.

Figures:
- `plot_envelope.png`: G(tau), nekStab against the reference curve.
- `plot_baseflow.png`: streamwise velocity of the base flow.
- `plot_optimal_perturbation.png`, `plot_optimal_response.png`: streamwise
  velocity of the optimal initial perturbation (`pRebfs0.f00001`) and of
  its response at t = tau (`orebfs0.f00001`), drawn from the tau = 57.905 run.
  The energy of the response divided by the energy of the perturbation is
  61688 (mass-weighted, GLL quadrature); G is 61690.

## Reference

- Blackburn, H. M., Barkley, D. & Sherwin, S. J. (2008). Convective
  instability and transient growth in flow over a backward-facing step.
  *J. Fluid Mech.* 603, 271-304.
