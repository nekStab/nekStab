# Validated Examples

The following examples are verified to run to completion and reproduce
the expected qualitative result with the canonical mesh and parameters
shipped in this repository.

## Cylinder (canonical 2128-element mesh)

| Case | Re | Mode | Expected result |
|---|---:|---|---|
| `baseflow/boostconv` | 100 | BoostConv (1.2) | converged steady state |
| `baseflow/newton` | 100 | Newton-GMRES (2.0) | converged steady state |
| `baseflow/newton_dyn` | 100 | Newton + `ifdyntol` (2.0) | converged steady state |
| `baseflow/sfd` | 50 | SFD baseline (1.1) | converged steady state |
| `baseflow/sfd_dyn` | 50 | SFD + dyn-tol (1.1) | converged steady state |
| `baseflow/sfd_dyn_oifs` | 50 | SFD + OIFS (1.1) | converged steady state |
| `stability/direct` | 100 | Direct LNSE (3.1) | leading eigenvalue near σ = 0.125 + 0.734i |
| `stability/adjoint` | 100 | Adjoint LNSE (3.2) | matches direct spectrum |
| `stability/animate_modes` | 100 | Mode animation (4.50) | snapshots per period |
| `postproc/sensitivity_budget_wavemaker` | 100 | Wavemaker (4.0) | wm, sr/si, pr/pi, tr/ti fields |
| `postproc/steady_force_sensitivity` | 100 | Steady force (4.41) | sensitivity field |
| `dns` | 100 | DNS (0.0) | vortex shedding |
| `ci_test` | 100 | DNS smoke (0.0) | 10-step smoke |
| `moving_cylinder` | 100 | DNS forced osc (0.0) | DNS subcase |

## Other geometries

| Case | Re / Ra | Mode | Expected result |
|---|---:|---|---|
| `lid_driven` | 3600 | Newton (2.0) | converged steady state |
| `poiseuille_OTD` | 5000 | OTD smoke | builds and runs |
| `thersyphon/baseflow` | Ra = 500 | Newton + thermal (2.0) | converged thermal steady state |
| `thersyphon/stability/direct` | Ra = 500 | Direct + thermal (3.1) | leading real eigenvalue (pitchfork) |
| `naca0012` | 2000 | Newton (2.0) | converged steady state |
| `flip_flop/baseflow` | 62 | Newton-UPO (2.1) | converged UPO |
| `flip_flop/stability/direct_Floquet` | 62 | Floquet direct (3.11) | leading multiplier just above unit circle |
| `tpjet/baseflow/newton` | 1900 | Forced UPO (2.2) | converged forced UPO |
| `tpjet/stability/direct_Floquet` | 1900 | Floquet direct (3.11) | leading multiplier above unit circle |

## Reproducing

Each validated case includes a `README.md` describing physics, mode
selection, and expected output, plus the necessary `.par`, `.usr`,
`SIZE`, and (where relevant) base-flow files.

```bash
cd example/<path>
mks <casename>
nekbmpi <casename> <N>     # or sbatch <case>.slurm on a cluster
python3 plot.py            # regenerate the evidence figure
```

## Deferred for v2.1

The following cases ship with their source and parameters but were
not part of the v2.0 verification sweep:

- `cylinder/baseflow/{newton_dyn_temp, newton_dyn_upo}`
- `cylinder/stability/{direct_Floquet, adjoint_Floquet, animate_modes_with_UPO}`
- `cylinder/{otd, modal, RANS}`
- `back_fstep/{baseflow, transient_growth}`
- `cubic_cavity`, `cubic_cavity_upo`
- `tpjet/baseflow/tdf`
- `slot_FST`
- `blasius`, `poiseuille_RANS`
