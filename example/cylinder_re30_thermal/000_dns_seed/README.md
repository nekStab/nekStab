# cylinder_re30_thermal / 000_dns_seed

DNS of the heated-cylinder wake at Re = 25, Ri = 0.2 (`userParam06`), started
from the uniform field and run for 100 convective times. The last snapshot was
the Newton seed of an earlier version of `../210_baseflow_newton/`.

The current Newton stage does not read this seed. It cold-starts from the
uniform field (`startFrom` is commented out in its `1cyl.par`), because the
earlier seed was written on a different mesh. This stage therefore does not
feed the chain; it shows how to produce a thermal seed.

Run: `mks 1cyl`, then `sbatch run.local.slurm`.
