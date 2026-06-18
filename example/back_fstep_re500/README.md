# Backward-Facing Step at Re=500

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This case is the backward-facing-step transient-growth example. The base flow
is used to study finite-time amplification in a stable but non-normal open
shear flow, so the important downstream stage is `330_transient_growth` rather
than a modal-instability threshold. Thesis source: chapter 2,
`Validation & Verification`, subsection `Backward-facing step`; AMR_Krylov_V5
Main `S5_examples.tex`, subsection `Backward-facing step`.

## Stages

- `210_baseflow_newton/` computes the steady fixed point used by stability
  analysis. Two-phase recipe (see its README): a short DNS settles the steady
  recirculation on the actual mesh (`bfs_dns.par`), then Newton-GMRES polishes
  that native seed to machine zero (`bfs_newton.par`). DNS seed generation is
  folded into this stage — there is no separate `000_dns/`.
- `330_transient_growth/` computes optimal perturbation and response over a
  finite horizon.

