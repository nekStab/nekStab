# cylinder_re100 — canonical full-pipeline showcase

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This is the canonical circular-cylinder wake ladder for fixed-point,
direct/adjoint, sensitivity, transient-growth, OTD, and snapshot-modal
workflows. The thesis uses the circular-cylinder wake to illustrate primary
instability, sensitivity, and secondary Floquet instability concepts. Thesis
source: chapter 2, `sec:ex:1cyl`, `sensitivity_cylinder`, and the following
pitchfork-bifurcation subsection.

**Target**: Re=100. Reynolds is fixed across all stages.

**Status**: scaffolded only — no case files yet. See per-stage READMEs.

**Pipeline** (16 stages):
000, 110 (x4 SFD variants), 120, 130, 140, 210 (x2 Newton variants), 310 (x2 direct variants), 311, 320, 321, 330, 411, 500, 600, 610, 620

**Cross-stage comparison** (the pedagogical payoff for 110): once 110 variants are wired, run
`scripts/nstab_compare.py example/cylinder_re100/110_baseflow_sfd/` to overlay residual decays.
