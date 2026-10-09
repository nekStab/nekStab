# flip_flop_re62 / 000_dns

**Stage**: DNS (userParam01=0), cold start.
**Re**: 62
**Ranks**: 8

Nonlinear reference trajectory / seed for downstream stages. Run with `sbatch run.local.slurm`.

## Reference
`ref/reference.json` checks the last snapshot, `2cyl0.f00010`: its time stamp (t = 100) and the norms of `vx` (1 %) and `vy` (2 %). The flow is still developing at t = 100, so the check is loose.
