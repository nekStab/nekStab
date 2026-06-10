# Cubic Cavity at Re=1914 — primary Hopf criticality (Newton-GMRES baseflow)

## Why Re=1914

The 3D cubic lid-driven cavity (aspect ratio Λ=1) undergoes its primary
**Andronov-Poincaré-Hopf bifurcation** at Re_c ≈ 1914 (subcritical;
Loiseau, Robinet & Leriche, *Fluid Dynamics Research* 48:6, 061421,
2016). The thesis numerically refines this value to **Re_c = 1916.63**
via Krylov-Schur on the time-stepper (Table 4.1 of Schuh-Frantz 2022,
chapter 4 §LinearStability).

Independent verifications:

- Feldman & Gelfgat (2010): Re_c = 1914
- Kuhlmann & Albensoeder (2014): Re_c ≈ 1919.5
- Liberzon (2011) experiment: 1700 ≤ Re_c ≤ 1970, St = 0.097
  (matches our St = 0.0934)

**At Re = 1914 the steady state still exists** but is on the verge of
instability — the leading eigenvalue pair sits on the imaginary axis
by construction. This is the canonical reference for stability
analysis: direct, adjoint, and structural sensitivity all live here.

## Physics

- 3D lid-driven cubic cavity, Λ=1 (cubic), Re=1914
- `viscosity = -1914.0` in `cav.par`
- `userParam01 = 2.0` — Newton-GMRES fixed point
- `userParam07 = 100` — Krylov subspace dimension

## Initial-condition workflow

The dir was renamed from `cubic_cavity_re2500/` on 2026-05-21 to align
with the Loiseau-Robinet-Leriche bifurcation reference. The shipped
`BF_cav0.f00001` is therefore the **Re = 2500 baseflow** from the prior
operating point and serves as a seed (Newton can converge unstable
steady states at any Re inside the basin):

1. Use `BF_cav0.f00001` (Re=2500) as the initial guess via
   `startFrom = BF_cav0.f00001` (already set in `cav.par`).
2. Run Newton-GMRES with `viscosity = -1914`. Expected: ~10 Newton
   iterations from this seed to converge the Re=1914 unstable steady
   state.
3. Save the resulting `BF_cav0.f00001` over the Re=2500 seed once
   convergence is verified.

## Status

`NEEDS_DNS_SEED` until the Newton-GMRES converges the Re=1914
baseflow. After convergence, downstream cases in this dir
(`310_stability_direct/`, `320_stability_adjoint/`,
`411_postproc_wavemaker/`, etc.) can run at the criticality reference.

## Relationship to cubic_cavity_re1950

- **`cubic_cavity_re1914/`** (this dir): primary Hopf bifurcation
  criticality reference. Used for stability analysis at Re_c.
- **`cubic_cavity_re1950/`** (sibling): just-supercritical UPO on the
  LC1 limit cycle (the limit cycle component of the post-Hopf
  intermittent attractor; T ≈ 10.77, St ≈ 0.093). Used for Floquet
  analysis on the periodic regime.

The two operating points together exercise both branches of the
post-bifurcation analysis: linear-stability AT criticality (Re=1914)
and Floquet of the periodic regime just past criticality (Re=1950).

## References

- Loiseau, J-Ch., Robinet, J-Ch., & Leriche, E. (2016). *Intermittency
  and transition to chaos in the cubical lid-driven cavity flow.*
  Fluid Dynamics Research 48(6), 061421.
- Schuh-Frantz thesis (2022), chapter 4 §`cavsubsec:LinearStability`,
  Table 4.1: Re_c = 1916.63, St = 0.0934 at Λ=1.0.
- Feldman, Y. & Gelfgat, A. (2010). *Oscillatory instability of a
  three-dimensional lid-driven flow in a cube.* Physics of Fluids 22,
  093602.
- Kuhlmann, H. C. & Albensoeder, S. (2014). *Three-dimensional flow
  instabilities in a lid-driven cavity.* Physics of Fluids 26, 042103.
- Liberzon, A., Feldman, Y., & Gelfgat, A. Yu. (2011). *Experimental
  observation of the steady-oscillatory transition in a cubic
  lid-driven cavity.* Physics of Fluids 23, 084106.
