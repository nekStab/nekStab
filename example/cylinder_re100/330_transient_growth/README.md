# Cylinder wake: optimal transient growth (uparam01 = 3.3)

## What this case shows

For a horizon tau, nekStab finds the initial perturbation q0 with the
largest energy gain

    G(tau) = max over q0 of  E(tau) / E(0)

as the leading eigenvalue of M^H M, where M = exp(tau L) is the
linearized flow map and M^H its adjoint. Each Krylov step is one forward
and one adjoint integration over tau.

The case runs at **Re = 40**, below the onset of vortex shedding
(Re_c = 46.6). Every eigenvalue of the steady wake decays, so any gain
above 1 is non-modal growth: a perturbation near the cylinder is
amplified while the flow carries it downstream (Cantwell & Barkley 2010).

## Setup

- Mesh: the shared cylinder mesh, 2128 elements, x in [-25, 100],
  y in [-25, 25].
- `viscosity = -40.0` (Re = 40). `endTime = 40` is the horizon tau.
- `userParam07 = 50` forces exactly one Krylov-Schur restart, so this
  stage also tests the restart path (see the comment in `1cyl.par`).
- Sponge `userParam08/09/10 = 5 / 5 / 1.7`: it damps the perturbation in
  the last 5 units before the outlet, in the direct and the adjoint run.
- `dt = 0.005` fixed; 8 MPI ranks.

## Base flow

`BF_1cyl0.f00001` is the steady Newton solution at Re = 40. Make it with
`../210_baseflow_newton/fp/` after you set `viscosity = -40.0` in its
`1cyl.par`: Newton converges in 5 iterations to 9.8e-10 (73 min on
8 ranks). Copy the result here. Field files are not tracked.

Run: `sbatch run.local.slurm`. G is the first value in
`Spectre_Hp_conv.dat`. The optimal perturbation is `pRe1cyl0.f00001`,
its response at t = tau is `ore1cyl0.f00001`.

## Result

One run per horizon (`growth_sweep.csv`, `plot_envelope.png`):

| tau | G(tau) |
|---:|---:|
| 10 | 18.0 |
| 20 | 75.4 |
| 40 | 275.0 |
| 60 | 384.7 |
| 80 | 345.8 |

The gain peaks between tau = 60 and 80 and then decays, as it must for a
stable flow. Wall time at tau = 40: 3 h 55 min on 8 ranks. These are
nekStab values on this mesh and domain; no published G(tau) at Re = 40
on the same domain was compared.

## Reference

- Cantwell, C. D. & Barkley, D. (2010). Computational study of subcritical
  response in flow past a circular cylinder. *Phys. Rev. E* 82, 026315.
