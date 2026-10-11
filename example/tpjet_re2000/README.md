# Harmonically forced jet (tpjet_re2000)

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

Axisymmetric jet at Re = 2000 (`viscosity = -2000`) with a harmonically forced inflow, `1 + 0.05 cos(2 pi St t)` at
St = 0.6, over a tanh shear layer of momentum thickness 0.025. The forced flow has a periodic orbit of period
T = 1 / St = 1.66667. The mesh has 4800 elements (160 by 30) and polynomial order 5 (`lx1 = 6`), the mesh of the
AMR_Krylov paper. In that paper the orbit becomes unstable at Re = 1371 through a period doubling (a real Floquet
multiplier through -1), and its Floquet spectrum at Re = 2000 has one real multiplier outside the unit circle.
The case is the one of Shaabani-Ardali, Sipp and Lesshafft (J. Fluid Mech. 2019). Sources: AMR_Krylov_V5/Main
`S5_examples.tex` section `sec:ex:hjet`; thesis chapter 6, `jet_sec:probform`, `jet_sec:global`, and
`jet_sec:secondary_jet`.

An orbit that is unstable at Re = 2000 cannot be found by marching the flow in time. The two methods that find it
are time-delayed feedback (TDF), which holds the flow on the orbit with a force proportional to the difference
between the flow and the flow one period earlier, and Newton-GMRES. Both stages start from the orbit at Re = 1900.

## Stages
| Stage | Mode | What it does | Result |
|---|---|---|---|
| `140_baseflow_tdf` | 1.4 | TDF from the orbit at Re = 1900; the flow approaches the orbit at Re = 2000 | residual 6.5e-8 at t = 600, 1.9e-3 from the Newton orbit |
| `210_baseflow_newton` | 2.2 | Newton-GMRES from the same seed; the orbit at Re = 2000 | 5 iterates, squared residual 3.5e-6 to 9.7e-9 |
| `311_stability_direct_floquet` | 3.11 | Floquet multipliers of that orbit | mu = -1.1577 (period doubling), sigma = 0.0878 |
| `321_stability_adjoint_floquet` | 3.21 | Floquet multipliers of the adjoint map of that orbit | mu = -1.1573; adjoint mode at the nozzle lip |

The seed `BF_Re1900_tpjet0.f00001` is the orbit at Re = 1900 of the AMR_Krylov paper, stored at polynomial order 5.
Stages 140 and 210 ship it. Stages 311 and 321 ship the orbit that stage 210 writes.

Each stage is self-contained (case files, start field, `ref/` baseline). See each stage's README and run
`scripts/check_against_ref.py <stage>` to verify a fresh run. All stages need at least 16 ranks (`lpmin = 16`).
