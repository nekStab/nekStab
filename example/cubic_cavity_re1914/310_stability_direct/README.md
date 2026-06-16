# cubic_cavity_re1914 / 310_stability_direct

**Stage**: 310_stability_direct
**uparam01**: 3.1
**Re**: 1914
**Start from**: BF_cav0.f00001
**Ranks**: 16

Stage rationale is shared in the repo-root `EXAMPLES.md`. Run with
`sbatch run.local.slurm`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.
