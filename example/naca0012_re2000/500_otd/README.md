# naca0012_re2000 / 500_otd

**Stage**: 500_otd
**uparam01**: 5
**Re**: 2000
**Start from**: cold start (useric)
**Ranks**: 8

Stage rationale is shared in the repo-root `EXAMPLES.md`. Run with
`sbatch run.local.slurm`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.
