# cubic_cavity_re1920 / 210_baseflow_newton/fp

**uparam01**: 2.0 (steady Newton-GMRES, k_dim = 100)
**Re**: 1920, 16 ranks (`lpmin = 16`)

Steady base flow of the cubic lid-driven cavity at Re = 1920, uniform lid. Re = 1920 lies above the
critical Reynolds number Re_c = 1916.63 of the thesis (Table 4.1). The steady flow exists and is
unstable, so a DNS cannot hold it and Newton-GMRES has to find it.

**Seed**: the last field of `../../000_dns/` (DNS at Re = 1900 from rest to t = 100), copied here as
`BF_cav0.f00001`. `startFrom = ... time=0` resets the time. Newton writes the converged flow over the
same file name. The seed is not steady (squared residual 8.5e-6, 85 000 times the tolerance of 1e-10).

**Result**: 12 iterations from a squared residual of 8.5e-6 to 9.7e-11 (`residu_newton.dat`), 122
Arnoldi vectors in total, 30 min on 40 ranks (Jean Zay, 2026-10-10). The residual is not monotone: it
rises once, from 3.5e-7 at iteration 5 to 8.5e-7 at iteration 6, and then falls by a factor of 2 to 7
per iteration. The converged `BF_cav0.f00001` is the base flow of `../../310_stability_direct/`.

**Reference**: `ref/reference.json` checks the first residual (8.5e-6, 5 %), the iteration count (12,
tolerance 2), the final residual (below the 1e-10 tolerance) and the norms of `vx`, `vy` and `vz` of the
base flow (1e-3). The first residual and the iteration count fail if the run starts from a converged
field. Verify with `scripts/check_against_ref.py example/cubic_cavity_re1920/210_baseflow_newton/fp`.

**Figure**: `python3 plot.py` writes `plot_newton.png`: (a) the residual history and (b) the speed of
the base flow on the plane z = 0. `ref/newton.png` is the figure of the run named in
`ref/reference.json`.

**Run**: `sbatch run.local.slurm`. Field files are not tracked.
