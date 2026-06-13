# lid_driven_re3600 / 000_dns

**Stage**: 000_dns
**uparam01**: 0
**Re**: 3600
**Start from**: cold start (useric)
**Ranks**: 4

Stage rationale is shared in the repo-root `EXAMPLES.md`. Run with
`sbatch run.local.slurm`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.
