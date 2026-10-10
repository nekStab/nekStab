# cubic_cavity_re1914 / 310_stability_direct

**uparam01**: 3.1 (Krylov-Schur, Arnoldi window endTime = 2.0, k_dim = 120,
6 eigenvalues converged)
**Re**: 1914, 16 ranks (`lpmin = 16`)

Copy the Newton `BF_cav0.f00001` of `../210_baseflow_newton/fp/` here, then
`sbatch run.local.slurm`. Converged eigenvalues: `Spectre_NSd_conv.dat` (sigma, omega in
units of U/L). The leading pair has sigma = -0.000183 and omega = 0.587 (St = 0.0934).
Result and comparison with the other stage: see `../README.md`.

**Run** (Jean Zay, 2026-10-10, 40 ranks): 2437 s. The six converged eigenvalues agree with those of a
16-rank run to 1e-8 in sigma and 1e-7 in omega.

**Reference**: `ref/reference.json` checks sigma and omega of the leading pair (tolerances 1e-3 and 2e-3).
Verify with `scripts/check_against_ref.py example/cubic_cavity_re1914/310_stability_direct`.

**Figure**: `python3 plot.py` writes `plot_spectrum.png`, the direct and the adjoint spectrum (it reads
`../320_stability_adjoint/Spectre_NSa_conv.dat`). `ref/spectrum.png` is the figure of that run.
