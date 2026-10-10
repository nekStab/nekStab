# cylinder_re100 — canonical full-pipeline showcase

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This is the canonical circular-cylinder wake ladder for fixed-point,
direct/adjoint, sensitivity, transient-growth, OTD, and snapshot-modal
workflows. The thesis uses the circular-cylinder wake to illustrate primary
instability, sensitivity, and secondary Floquet instability concepts. Thesis
source: chapter 2, `sec:ex:1cyl`, `sensitivity_cylinder`, and the following
pitchfork-bifurcation subsection.

**Target**: Re = 100. Two stages differ: `120_baseflow_boostconv` runs at Re = 50 (the regime the
BoostConv method was designed for) and `330_transient_growth` at Re = 40 (below the onset of shedding,
Re_c = 46.6). The Reynolds number of a stage is the (negative) `viscosity` in its `1cyl.par`.

**Stages** (each folder has its own README):

| stage | what it computes |
|---|---|
| `000_dns/` | DNS of the shedding wake; `ci_test/` is the 10-step smoke test of the CI |
| `110_baseflow_sfd/` | steady base flow by selective frequency damping, six variants: the Akervik and the Casacuberta form, each plain, with dynamic tolerance (`_dyn`), and with dynamic tolerance and OIFS time stepping (`_dyn_oifs`) |
| `120_baseflow_boostconv/` | steady base flow by BoostConv |
| `130_baseflow_dmt/` | steady base flow by dynamic mode tracking |
| `210_baseflow_newton/` | steady base flow by Newton-GMRES, `fp` (fixed tolerance) and `dyn` (dynamic tolerance) |
| `140_baseflow_meanflow/` | mean flow of the shedding wake, averaged over 25 whole periods; the DNS stops when the average closes |
| `310_stability_direct/direct/` | direct eigenvalue problem of the steady wake |
| `315_stability_direct_meanflow/` | direct eigenvalue problem of the mean flow: growth rate near zero, frequency of the shedding |
| `320_stability_adjoint/` | adjoint eigenvalue problem |
| `330_transient_growth/` | optimal transient growth of the stable wake at Re = 40 |
| `410_postproc_animate_modes/`, `411_postproc_wavemaker/` | mode animation, wavemaker |
| `500_otd/` | optimally time-dependent modes about the steady wake |
| `600_modal_pod/`, `610_modal_dmd/`, `620_modal_spod/` | POD, DMD and SPOD of the shedding wake |

**Cross-stage comparison**: `scripts/nstab_compare.py example/cylinder_re100/110_baseflow_sfd/` overlays
the residual decays of the SFD variants.
