# Triple-port jet (tpjet)

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This axisymmetric forced-jet case supports periodic-orbit and Floquet analysis
of vortex-pairing dynamics. The AMR_Krylov examples discuss the harmonically
forced jet as a period-doubling Floquet case; the named thesis jet chapter
covers jet transition, nonlinear regimes, base-flow description, global
stability, and secondary instability. Sources: AMR_Krylov_V5/Main
`S5_examples.tex` section `sec:ex:hjet`; thesis chapter 6,
`jet_sec:probform`, `jet_sec:global`, and `jet_sec:secondary_jet`.

Axisymmetric triple-port jet — forced-orbit base flow and Floquet stability.

## ⚠ Directory naming vs operating point

This directory is named `tpjet_re2005`, but the operating points of the shipped
stages are **mixed**:

| Stage | Mode | Re |
|---|---|---|
| `210_baseflow_newton/forced_po` | 2.2 forced PO | **1900** |
| `311_stability_direct_floquet`  | 3.11 Floquet  | **1900** |
| `baseflow/tdf` (not yet migrated) | 1.4 TDF | 2005 |
| `130_baseflow_dmt` (config only) | 1.3 DMT | 1900 |

The migrated, data-complete stages are at **Re = 1900** (above Re_c ≈ 1371). The
only Re = 2005 stage (`tdf`) has no converged result yet, so it is not migrated.
**Recommendation:** rename this directory to `tpjet_re1900` (or split the TDF
case into its own `tpjet_re2005`) so the name matches the operating point.

## Stages
- `210_baseflow_newton/forced_po` — forced periodic orbit via Newton-GMRES.
- `311_stability_direct_floquet` — direct Floquet stability of that orbit.

Each stage is self-contained (config + IC + `ref/` baseline). See each stage's
README and run `scripts/check_against_ref.py <stage>` to verify a fresh run.
