# cubic_cavity_re1914 / 330_transient_growth

**uparam01**: 3.3 (optimal transient growth), Re = 1914, 16 ranks

Optimal transient growth of the steady cubic-cavity flow at Re = 1914. The flow is linearly
stable (leading sigma = -1.8e-4, `../310_stability_direct/`), yet a well-chosen initial
perturbation gains a factor 229 in energy over the horizon tau = 10 (units L/U). For a horizon tau,
nekStab finds the largest energy gain

    G(tau) = max over q0 of  E(tau) / E(0)

as the leading eigenvalue of M^H M, where M = exp(tau L) is the linearised flow map and M^H its
adjoint. Each Krylov step is one forward and one adjoint integration over tau.

## Setup

- Base flow: `BF_cav0.f00001`, the Newton solution of `../210_baseflow_newton/fp/`. Copy it here;
  field files are not tracked.
- `endTime = 10.0` is the horizon tau. `userParam07 = 30` is the size of the Krylov basis: one
  Arnoldi pass of 30 vectors converges 13 eigenvalues of M^H M, and `userParam03 = 2` writes the
  two leading optimal perturbations.
- Isothermal, no sponge, uniform lid as in the thesis. `dt` follows from `targetCFL = 0.5` and is
  then adjusted so that the forward and the adjoint integration cover exactly tau (6958 steps of
  1.437e-3 for tau = 10).
- 16 MPI ranks (`lpmin = 16`). One Arnoldi iteration takes 3.7 min on 40 ranks of Jean Zay, so
  the run takes 1 h 50 min there. `run.local.slurm` asks for 8 h.

Run: `sbatch run.local.slurm`. G is the first value in `Spectre_Hp_conv.dat`. The optimal
perturbation is `pRecav0.f00001` and its response at t = tau is `orecav0.f00001`.
Plot: `python3 plot.py`.

## Result

G(10) = 229.19 and the second gain is 215.96. The stored perturbation has unit energy. The
response has the energy 228.50, so the gain of the stored perturbation is 0.3 % below G. G is the
Rayleigh quotient of M^H M, which uses the adjoint integration; the energy ratio of the fields uses
the forward integration only. The two agree exactly only when the adjoint integration is the
transpose of the forward one, so the 0.3 % measures how far it departs from it on this mesh.

On the plane z = 0 the optimal perturbation is largest at (x, y) = (-0.02, 0.41), in the shear
layer under the lid. At t = 10 the response is largest at (-0.41, -0.24), near the left wall, on
the far side of the main vortex. No published G(tau) of the cubic cavity at this Reynolds number
was compared.

`growth_sweep.csv` holds G for the horizons that were run as copies of this folder with another
`endTime`. `python3 plot.py` draws G(tau) with the sweep (a) and the speed of the optimal
perturbation (b) and of its response (c) on the plane z = 0.

## Reference
`ref/reference.json` checks the two leading gains (1e-4), the energy of the stored perturbation
(1 within 1e-5) and the gain of the stored perturbation divided by G (within 0.005 of 0.997).
The last two come from `growth_check.dat`, which `plot.py` writes from the field files. The last
check fails when the stored response is not exp(tau L) applied to the stored perturbation: a
response scaled by G gives the ratio G, not 1.

## Figure
`ref/growth.png` is the figure of the run named in `ref/reference.json`.
