# lid_driven_re3600 / 320_stability_adjoint

**Stage**: 320_stability_adjoint · **uparam01**: 3.2 · **Re**: 3600 · **Ranks**: 4
**Start from**: `BF_cav0.f00001`, the Newton base flow of `../210_baseflow_newton/fp` (copy the
result `BF_cav0.f00001` here). Stage rationale: repo-root `EXAMPLES.md`.

Adjoint eigenvalue problem (Krylov-Schur, `k_dim = 120`) about the same base flow as
`../310_stability_direct/`, `lx1 = 12`. The adjoint operator has the spectrum of the direct one,
so this stage checks the adjoint: at `lx1 = 12` the leading eigenvalue, −0.016796, is within 4e-6
of the direct −0.016793, and the pair is −0.04197 ± 0.77165 i against −0.04198 ± 0.77165 i.
`ref/reference.json` compares the adjoint eigenvalues with the direct ones.
