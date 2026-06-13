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

- `000_dns/` stores the nonlinear DNS scaffold.
- `210_baseflow_newton/` computes the fixed point used by stability analysis.
- `330_transient_growth/` computes optimal perturbation and response over a
  finite horizon.

