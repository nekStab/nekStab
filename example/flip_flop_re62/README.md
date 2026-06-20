# Side-by-Side Cylinders at Re=62

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This is the side-by-side-cylinder flip-flop case. It is organized around a
periodic orbit and Floquet analysis of the secondary Neimark-Sacker mechanism
associated with the flip-flop dynamics. Thesis source: chapter 2,
`sec:ex:2cyl` and `Neimark-Sacker bifurcation`; AMR_Krylov_V5 Main
`S5_examples.tex`, subsection `The flow past side-by-side circular cylinders`.

## Stages

- `000_dns/` supplies the nonlinear trajectory.
- `210_baseflow_newton/` computes the periodic orbit.
- `311_stability_direct_floquet/` and `321_stability_adjoint_floquet/` analyze
  the Floquet instability.
- `411_postproc_wavemaker/` localizes the coupled direct/adjoint sensitivity.

