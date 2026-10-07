# cubic_cavity_re1914 / 210_baseflow_newton/fp

**uparam01**: 2.0 (steady Newton-GMRES, k_dim = 100)
**Re**: 1914, 16 ranks

Steady base flow of the cubic lid-driven cavity at Re = 1914, uniform lid.

**Seed**: the last field of `../../000_dns/` (DNS from rest to t = 100), copied
here as `BF_cav0.f00001`. `startFrom = ... time=0` resets the time. Newton
writes the converged flow over the same file name.

**Result**: 12 iterations to a residual of 8.4e-11 (`residu_newton.dat`),
76 min on 16 ranks. The converged `BF_cav0.f00001` is the base flow of
`../../310_stability_direct/` and `../../320_stability_adjoint/`.

**Run**: `sbatch run.local.slurm`. Field files are not tracked.
