# lid_driven_re3600 / 310_stability_direct/direct

**Stage**: 310_stability_direct (direct)
**uparam01**: 3.1
**Re**: 3600
**Start from**: BF_cav0.f00001
**Ranks**: 4

Stage rationale is shared in the repo-root `EXAMPLES.md`. Run with
`sbatch run.local.slurm`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.

**KNOWN ISSUE (2026-06-13)**: stability run aborts at "Mesh check failed" — the
lid_driven mesh/SIZE that works for the root baseflow fails the stability stage
check. Same class as thermal dyn_temp. Parked; needs mesh/SIZE alignment.
