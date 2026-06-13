# lid_driven_re3600 / 330_transient_growth

**Stage**: 330_transient_growth
**uparam01**: 3.3
**Re**: 3600
**Start from**: BF_cav0.f00001
**Ranks**: 4

Stage rationale is shared in the repo-root `EXAMPLES.md`. Run with
`sbatch run.local.slurm`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.
