# cubic_cavity_re1914 / 320_stability_adjoint

**uparam01**: 3.2 (Krylov-Schur, Arnoldi window endTime = 2.0, k_dim = 120,
6 eigenvalues converged)
**Re**: 1914, 16 ranks (`lpmin = 16`), 25-40 min

Copy the Newton `BF_cav0.f00001` of `../210_baseflow_newton/fp/` here, then
`sbatch run.local.slurm`. Converged eigenvalues: `Spectre_NSa_conv.dat` (sigma, omega in
units of U/L). Result and comparison with the other stage: see `../README.md`.
The flow is stable at Re = 1914 (leading sigma = -0.0170).
