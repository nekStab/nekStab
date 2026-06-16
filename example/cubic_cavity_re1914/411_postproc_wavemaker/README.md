# cubic_cavity_re1914 / 411_postproc_wavemaker

**Stage**: 411_postproc_wavemaker
**uparam01**: 4.11
**Re**: 1914
**Start from**: BF_cav0.f00001
**Ranks**: 16

Stage rationale is shared in the repo-root `EXAMPLES.md`. Run with
`sbatch run.local.slurm`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.

**PARKED (Re=1914 steady)**: 411_postproc_wavemaker needs a periodic orbit / snapshot series that does not exist below the Hopf. N/A here; applies at re1950 (supercritical).
