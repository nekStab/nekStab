# cubic_cavity_re1914 / 320_stability_adjoint

**uparam01**: 3.2 (Krylov-Schur, Arnoldi window endTime = 2.0, k_dim = 120,
6 eigenvalues converged)
**Re**: 1914, 16 ranks (`lpmin = 16`)

Copy the Newton `BF_cav0.f00001` of `../210_baseflow_newton/fp/` here, then
`sbatch run.local.slurm`. Converged eigenvalues: `Spectre_NSa_conv.dat` (sigma, omega in
units of U/L). The leading pair has sigma = -0.000336 and omega = 0.587 (St = 0.0934).
Result and comparison with the other stage: see `../README.md`.
