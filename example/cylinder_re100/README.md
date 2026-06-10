# cylinder_re100 — canonical full-pipeline showcase

**Target**: Re=100. Reynolds is fixed across all stages.

**Status**: scaffolded only — no case files yet. See per-stage READMEs.

**Pipeline** (16 stages):
000, 110 (x4 SFD variants), 120, 130, 140, 210 (x2 Newton variants), 310 (x2 direct variants), 311, 320, 321, 330, 411, 500, 600, 610, 620

**Cross-stage comparison** (the pedagogical payoff for 110): once 110 variants are wired, run
`scripts/nstab_compare.py example/cylinder_re100/110_baseflow_sfd/` to overlay residual decays.

**Migration**: existing cylinder cases live at example/cylinder/. They will be linked or copied into the new stage directories in a later migration step.
