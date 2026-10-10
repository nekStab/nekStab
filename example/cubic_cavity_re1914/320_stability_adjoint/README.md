# cubic_cavity_re1914 / 320_stability_adjoint

**uparam01**: 3.2 (Krylov-Schur, Arnoldi window endTime = 2.0, k_dim = 120,
6 eigenvalues converged)
**Re**: 1914, 16 ranks (`lpmin = 16`)

Copy the Newton `BF_cav0.f00001` of `../210_baseflow_newton/fp/` here, then
`sbatch run.local.slurm`. Converged eigenvalues: `Spectre_NSa_conv.dat` (sigma, omega in
units of U/L). The leading pair has sigma = -0.000336 and omega = 0.587 (St = 0.0934).
Result and comparison with the other stage: see `../README.md`.

**Gate run** (Jean Zay, 2026-10-10, 40 ranks): 2534 s. The six converged eigenvalues agree with those of a
16-rank run to 1e-8 in sigma and 1e-7 in omega. The leading adjoint growth rate, -3.4e-4, differs from the
direct one, -1.8e-4, by 1.5e-4; the frequencies differ by 4e-4.

**Reference**: `ref/reference.json` checks sigma and omega of the leading pair (tolerances 1e-3 and 2e-3).
Verify with `scripts/check_against_ref.py example/cubic_cavity_re1914/320_stability_adjoint`.

**Figure**: `python3 plot.py` writes `plot_adjoint.png`: (a) the adjoint spectrum and (b) the amplitude of
the leading adjoint mode on the plane z = 0. `ref/adjoint.png` is the figure of the gate run.
