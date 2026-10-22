# lid_driven_re3600 / 310_stability_direct

**Stage**: 310_stability_direct · **uparam01**: 3.1 · **Re**: 3600 · **Ranks**: 4
**Start from**: `BF_cav0.f00001`, the Newton base flow of `../210_baseflow_newton/fp` (copy the
result `BF_cav0.f00001` here). Stage rationale: repo-root `EXAMPLES.md`.

Direct eigenvalue problem (Krylov-Schur, at most `k_dim = 120` vectors; the run stops when the 3 leading
values have converged, at step 114) of the regularized lid-driven cavity at
Re = 3600, `lx1 = 12`. The flow is linearly stable: the leading eigenvalue is real,
σ = −0.01679, and the next two form a complex pair −0.0420 ± 0.7717 i. `../320_stability_adjoint/`
gives the same spectrum; the regularized lid has no corner singularity, so the adjoint problem
is valid. `ref/reference.json` holds the values of this run. The spectrum depends on the
resolution (`../README.md`).

Run of 2026-10-09 (Jean Zay, `ifort`, 4 ranks, base flow of the 210 run of that day): wall time 619 s.
`plot.py` writes `plot_bf.png` (velocity magnitude of the base flow), `plot_spectrum.png` (converged eigenvalues,
growth rate against frequency) and `plot_mode.png` (vy of the leading direct mode); `ref/` holds copies as
`bf.png`, `spectrum.png` and `mode.png`.
