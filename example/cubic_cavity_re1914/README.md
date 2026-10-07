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

The case runs the 3D lid-driven cubic cavity at Re = 1914, just below the
primary Hopf bifurcation (Re_c = 1916.63 in the thesis, 1914-1919.5 in the
literature, St = 0.0934). The steady flow is stable by a small margin. The
leading eigenvalues (sigma, omega in units of U/L):

| mode | direct | adjoint |
|---|---|---|
| pair 1 (T1) | -0.000183 +- 0.58684i | -0.000336 +- 0.58723i |
| pair 2 | -0.02049 +- 0.56094i | -0.02080 +- 0.56108i |
| pair 3 | -0.02607 +- 0.14771i | -0.02598 +- 0.14838i |

Pair 1 has St = omega / 2 pi = 0.0934, the value of the thesis, and sits
0.0002-0.0003 below the imaginary axis. Direct and adjoint agree to 2e-4 in
sigma and 7e-4 in omega.

## Lid

The lid moves in x with u = 1 on the whole lid, as in the thesis: no
regularisation removes the singularity at the lid edges. Between 2026-06-16
and 2026-10 the cases used a regularized lid; the flow under it is more
strongly stable and the T1 pair is absent (leading sigma = -0.060 for
(1-(2x)^2)^2 (1-(2z)^2)^2, -0.017 for exponent 18). Do not regularize the
lid to reproduce the thesis.

## Workflow

1. `000_dns/`: DNS from rest at Re = 1914. The flow reaches a steady state
   by t = 100 (27 min on 16 ranks, 69 571 steps).
2. `210_baseflow_newton/fp/`: copy the last field of the DNS to
   `BF_cav0.f00001` and run Newton. It converges in 12 iterations to
   8.4e-11 (76 min on 16 ranks). `time=0` in `startFrom` resets the time.
   Starting Newton from a field of a different lid fails: the fixed time
   step taken from that field is too large (CFL 9).
3. `310_stability_direct/` and `320_stability_adjoint/`: copy the Newton
   `BF_cav0.f00001` into each folder and run. The Arnoldi window is
   endTime = 2.0 and k_dim = 120. About 95 min (direct) and 146 min
   (adjoint) on 16 ranks. `python 310_stability_direct/plot.py` draws both
   spectra.

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
