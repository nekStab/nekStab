# cubic_cavity_re1914 / 210_baseflow_newton/fp

**uparam01**: 2.0 (steady Newton-GMRES, k_dim = 100)
**Re**: 1914, 16 ranks

Steady base flow of the cubic lid-driven cavity at Re = 1914, uniform lid.

**Seed**: the last field of `../../000_dns/` (DNS from rest to t = 100), copied
here as `BF_cav0.f00001`. `startFrom = ... time=0` resets the time. Newton
writes the converged flow over the same file name.

**Result**: 12 iterations from a squared residual of 9.0e-6 to 8.4e-11 (`residu_newton.dat`),
59 min on 16 ranks (Jean Zay, gate run 2026-10-09). The converged `BF_cav0.f00001` is the base flow of
`../../310_stability_direct/` and `../../320_stability_adjoint/`.

**Reference**: `ref/reference.json` checks the first residual (9.0e-6, 5 %), the iteration count (12,
tolerance 2) and the final residual (below the 1e-10 tolerance). The first residual fails if the run
starts from a converged field.

**Figure**: `python3 plot.py` writes `plot_newton.png`: (a) the residual history and (b) the speed of
the base flow on the plane z = 0. `ref/newton.png` is the figure of the gate run.

**Run**: `sbatch run.local.slurm`. Field files are not tracked.
