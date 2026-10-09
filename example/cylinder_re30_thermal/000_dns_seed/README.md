# cylinder_re30_thermal / 000_dns_seed

DNS of the heated-cylinder wake at Re = 25, Ri = 0.2 (`userParam06`), started
from the uniform field and run for 100 convective times. The last snapshot was
the Newton seed of an earlier version of `../210_baseflow_newton/`.

The current Newton stage does not read this seed. It cold-starts from the
uniform field (`startFrom` is commented out in its `1cyl.par`), because the
earlier seed was written on a different mesh. This stage therefore does not
feed the chain; it shows how to produce a thermal seed.

Run: `mks 1cyl`, then `sbatch run.local.slurm` (8 ranks; 8 min on Jean Zay with `ifort`, 2026-10-09).
`plot.py` writes `plot_snapshot.png` (temperature of the last snapshot); `ref/` holds a copy as `snapshot.png`.

## Reference
`ref/reference.json` checks the last snapshot, `1cyl0.f00004`: its time stamp (t = 100) and the norms of `vx` and of the temperature (1 % each). The stage feeds no other stage, so the check is loose.
