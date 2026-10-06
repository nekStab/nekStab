# Cubic Cavity at Re=1914 — steady base flow and its stability

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This case represents the cubic lid-driven cavity near the primary oscillatory
instability of the steady 3D base flow. It is the steady-state counterpart to
`cubic_cavity_re1950`: compute the base flow, then use direct, adjoint,
transient-growth, wavemaker, modal, or OTD stages to interrogate the linear
mechanisms around that state. Thesis source: chapter 4,
`cav:sec:problem_formulation`, `cavsubsec:BaseFlows`, and
`cavsubsec:LinearStability`.

## Result in this repository

The case runs the 3D lid-driven cubic cavity at Re = 1914, the value of the
primary Hopf bifurcation in the literature (Re_c = 1914-1919.5, St = 0.0934,
omega = 0.587). **With the lid and mesh of this case the steady flow is
stable at Re = 1914.** The leading eigenvalues (sigma, omega in units of U/L)
are:

| mode | direct | adjoint |
|---|---|---|
| real | -0.01700 | -0.01708 |
| pair 1 | -0.03308 +- 0.4560i | -0.03306 +- 0.4560i |
| pair 2 | -0.05277 +- 0.4733i | -0.05295 +- 0.4735i |

Direct and adjoint agree to 5e-4 in sigma. The neutral pair near
omega = 0.587 does not appear, so this case does not reproduce the
published critical Reynolds number. Checked and ruled out: the mesh. At
polynomial order 7 (10^3 elements, 8^3 points each) the same modes are
sigma = -0.01678, -0.03284 +- 0.4562i, -0.05287 +- 0.4735i, within 2 %.
Not checked: the lid profile (see below), the box size convention, the
time-stepper filter. The case is a working example of 3D steady Newton plus
direct and adjoint stability, not a reproduction of the critical point.

## Lid

The lid moves in x with u = (1-(2x)^18)^2 (1-(2z)^18)^2 (Leriche & Gavrilakis,
2000). It is 1 except within about 5 % of the edges and zero at the edges, so
the adjoint problem has no corner singularity. Before 2026-10 the cases used
(1-(2x)^2)^2 (1-(2z)^2)^2, which is much weaker over the whole lid; the flow
under that profile was strongly stable (sigma = -0.060 for the leading mode).

## Workflow

1. `000_dns/`: DNS from rest at Re = 1914. The flow reaches a steady state
   by t = 100 (14 min on 16 ranks).
2. `210_baseflow_newton/fp/`: copy the last field of the DNS to
   `BF_cav0.f00001` and run Newton. It converges in 5 iterations to
   8e-11 (8 min on 16 ranks). `time=0` in `startFrom` resets the time.
   Starting Newton from a field of another lid fails: the fixed time step
   taken from that field is too large for the new flow (CFL 9).
3. `310_stability_direct/` and `320_stability_adjoint/`: copy the Newton
   `BF_cav0.f00001` into each folder and run. The Arnoldi window is
   endTime = 2.0 and k_dim = 120. Each takes about 25-40 min on 16 ranks.
   `python 310_stability_direct/plot.py` draws both spectra.

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
