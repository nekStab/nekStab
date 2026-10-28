# lid_driven_re3600 / 320_stability_adjoint

**Stage**: 320_stability_adjoint · **uparam01**: 3.2 · **Re**: 3600 · **Ranks**: 4
**Start from**: `BF_cav0.f00001`, the Newton base flow of `../210_baseflow_newton/fp` (copy the
result `BF_cav0.f00001` here). Stage rationale: repo-root `EXAMPLES.md`.

Adjoint eigenvalue problem (Krylov-Schur, at most `k_dim = 120` vectors; the run stops when the 3 leading
values have converged, at step 116) about the same base flow as
`../310_stability_direct/`, `lx1 = 12`. The adjoint operator has the spectrum of the direct one,
so this stage checks the adjoint: at `lx1 = 12` the leading eigenvalue, −0.016796, is within 4e-6
of the direct −0.016792, and the pair is −0.04191 ± 0.77156 i against −0.04192 ± 0.77156 i (base flow of the t = 100 seed, 2026-10-28).
`ref/reference.json` compares the adjoint eigenvalues with the direct ones.

Run of 2026-10-09 (Jean Zay, `ifort`, 4 ranks, base flow of the 210 run of that day): wall time 609 s.
`plot.py` writes `plot_bf.png` (velocity magnitude of the base flow), `plot_spectrum.png` (converged eigenvalues,
growth rate against frequency) and `plot_mode.png` (vx (a) and vy (b) of the leading adjoint mode, sign fixed so that vy is positive where |vy| is largest); `ref/` holds copies as
`bf.png`, `spectrum.png` and `mode.png`. The adjoint and direct modes look alike but differ: after matching norm and sign, the relative L2 difference is 0.30 in vx and 0.44 in vy (correlation 0.95 and 0.91).
