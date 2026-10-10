# Cubic Cavity at Re=1920 — unstable steady state

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This case runs the cubic lid-driven cavity just above the primary Hopf bifurcation
(Re_c = 1916.63 in the thesis, Table 4.1). The steady flow still exists and is unstable. A DNS
cannot hold it, so Newton-GMRES finds it, and a direct stability stage computes the leading
eigenpair about it. It is the supercritical neighbour of `cubic_cavity_re1914` (the steady flow
just below Re_c) and of `cubic_cavity_re1950` (the saturated limit cycle). Thesis source:
chapter 4, `cav:sec:problem_formulation`, `cavsubsec:BaseFlows` and `cavsubsec:LinearStability`.

## Lid

The lid moves in x with u = 1 on the whole lid, as in the thesis and in `cubic_cavity_re1914`.
The lid is not regularized.

## Workflow

1. `000_dns/`: DNS from rest at Re = 1900 to t = 100 (21 min, 69 571 steps on 16 ranks). Re = 1900
   lies below Re_c, so the flow settles towards its steady state. The last field is close to the
   steady state but not on it.
2. `210_baseflow_newton/fp/`: copy the last field of the DNS to `BF_cav0.f00001` and run Newton at
   Re = 1920. It converges in 12 iterations to 9.7e-11 (30 min on 40 ranks, Jean Zay).
   `time=0` in `startFrom` resets the time.
3. `310_stability_direct/`: copy the Newton `BF_cav0.f00001` here and run. The Arnoldi window is
   `endTime` = 2.0 and `k_dim` = 120.

## References

- Schuh-Frantz thesis (2022), chapter 4, Table 4.1: Re_c = 1916.63, St = 0.0934 at aspect
  ratio 1.
- Loiseau, J-Ch., Robinet, J-Ch., & Leriche, E. (2016). *Intermittency and transition to chaos in
  the cubical lid-driven cavity flow.* Fluid Dynamics Research 48(6), 061421.
