# Re=180 DNS: seed and period for the periodic orbit

A plain DNS (`userParam01 = 0`) at Re=180, started from rest. It gives the
two inputs of the Newton periodic-orbit stage `../210_baseflow_newton/`:

1. **Seed field.** The wake reaches its limit cycle by t = 150. The last
   checkpoint, `1cyl0.f00003` at t = 300, is copied to
   `../210_baseflow_newton/rstcyl0.f00001`.
2. **Period guess.** `1cyl.his` places one probe at (x, y) = (5, 0). Nek
   appends `time u v p` at each step. `python fft_period.py` takes the
   transverse velocity on t = 150-300 and gives the period two ways:
   FFT peak T = 5.19012, and a count of 28 zero crossings T = 5.19022.
   Newton converges to T = 5.189626 (St = 0.19269).

Re=180 is below the 3D mode A threshold (Re = 188.5). In 2D the wake has
no secondary instability, so the DNS settles on the same limit cycle that
Newton computes.

Run: `sbatch run.local.slurm` (8 ranks; 34 763 steps, 17 min on x2, a 6-core machine).
