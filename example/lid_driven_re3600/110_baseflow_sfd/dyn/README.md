# lid_driven_re3600 / 110_baseflow_sfd/dyn

**Stage**: 110_baseflow_sfd (dyn)
**uparam01**: 1.1
**Re**: 3600
**Start from**: BF_cav0.f00001
**Ranks**: 4

Stage rationale is shared in the repo-root `EXAMPLES.md`. Run with
`sbatch run.local.slurm`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.
