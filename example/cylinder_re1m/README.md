# Cylinder at Re=1e6

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This directory is the high-Reynolds-number cylinder. The point of the case is the gap between a base flow and a mean flow. That gap is small near onset and large far above it. Thesis source for the cylinder wake: chapter 2, `sec:ex:1cyl`. Method source: chapter 2, `Large-scale eigensolvers`.

## Two states

Stability is computed on a base flow and on a mean flow. Each state is either steady or periodic. One spectrum does not stand in for the other.

A base flow is a solution of the instantaneous equations. The steady base is a fixed point, from SFD or Newton. The periodic base is a periodic orbit. Stability of a steady base is the spectrum of the Jacobian at that field (`userParam01 = 3.1`). Stability of a periodic base is Floquet over one period (`userParam01 = 3.11`).

A mean flow is an average of a saturated unsteady run. The run may be DNS, filtered LES, or URANS. The steady mean is a time average. The periodic mean is a phase average. Stability is the spectrum about that average.

## Reynolds stress

The Reynolds stress is the fluctuation correlation

    R_ij = <u_i u_j> - <u_i><u_j>

It is the same object in DNS, filtered LES, and URANS. It is not the RANS closure. k, tau, and eddy viscosity are terms inside a URANS instantaneous equation. They are not R_ij. A URANS run can carry both.

On an exact fixed point the fluctuation is zero, so R_ij is zero. On an SFD trajectory that has not reached a fixed point, R_ij is the correlation of the residual SFD fluctuation. That number belongs to the base-flow example.

On the mean-flow example, R_ij is the time average, or the phase average, of the fluctuation products of the saturated run. That stress is the nonlinear distortion of the mean away from the base. It is small near the bifurcation and large far above it. Re=1e6 is the far case, which is why both examples are required.

## Reynolds-stress force

The force is off unless `nekStab_usrchk` calls `reynolds_enable`. That call is not a mode, so it works with DNS, SFD, Newton, stability, and Floquet. `nekStab_init` then loads `FRS<session>0.f00001` once and keeps that field for the run. `userf` adds it beside the other forces while `jp == 0`. A missing file logs and leaves the force off. The write example is `120_urans/1cyl.usr`: at the last step it calls `reynolds_commit`, which writes `RS1`, `RS2`, and `FRS` and does not arm them. One commit per run writes `FRS<session>0.f00001`. `310_stability_direct` does not call `reynolds_enable`.

## Examples

### Base, steady

`110_baseflow_sfd/` holds the SFD field `base_converged.f00001`. It is not a fixed point. The SFD residual floors near 5e-4, and transverse velocity is not zero. Figures: `110_baseflow_sfd/base_fields.png`, `110_baseflow_sfd/sfd_residual.png`.

`310_stability_direct/coupled/` and `quasilaminar/` linearize about that field. Those spectra are base-flow spectra of an imperfect steady base. Coupled versus quasilaminar is a choice about the closure perturbation. It is not the base-versus-mean choice, and it is not a Reynolds stress.

Explicit-filter results (`filterWeight = 1`, `filterCutoffRatio = 0.67`): coupled leads at σ = 3.41, St = 0.491, on the cylinder. Quasilaminar leads at σ = 1.15, St = 0, also on the cylinder. A wake-frequency mode remains in the quasilaminar spectrum at St = 0.189, σ = 0.261, and is not leading. The hpfrt result, σ = 0.237 at St = 0.194, is in `spectre_hpfrt.tar`.

Newton in `210_baseflow_newton/` did not produce a fixed point. Two iterations left the residual at 0.54, and the outposted velocity matches the SFD field. Figure: `210_baseflow_newton/newton_vs_sfd.png`.

Still required: a converged fixed point, then a new direct spectrum of that fixed point. The numbers above are not that spectrum.

### Base, periodic

Not started. Requires a periodic orbit of the instantaneous equations, then Floquet.

### Mean, steady

Not started. Requires a saturated unsteady run (DNS, filtered LES, or a URANS that fluctuates), the time average of that run as the linearization field, and R_ij from the same run.

`120_urans/` is not this example. The unforced k-tau march on this mesh does not shed. It leaves the SFD field and sits on a steady wake (`steady_wake.f00001`, mid-wake u about 0.61). Over a 20-time-unit window the resolved fluctuation correlation is about 10^{-3}. That is residual drift of a steady state, not the Reynolds stress of a saturated unsteady wake. A June run with v about ±0.87 was a different mesh. That figure is not kept.

### Mean, periodic

Not started. Requires a saturated periodic run, a phase average, and the phase-conditioned stress.

## Stages

- `000_dns/` is the nonlinear k-tau march on the 2D mesh (`lelg=1480`). Its `.par` is still the Re=40000 continuation point.
- `110_baseflow_sfd/` builds the SFD base the current stability stages restart from.
- `120_urans/` is the bare forward march. It documents that this operator does not fluctuate.
- `210_baseflow_newton/` is a failed Newton probe from the SFD field.
- `310_stability_direct/coupled/` is the coupled closure operator on the imperfect SFD base.
- `310_stability_direct/quasilaminar/` is the frozen-eddy-viscosity operator on the same base.

## Mesh

The run mesh is each stage's `1cyl.re2` (1480 elements, body-fitted, closest node at r = 0.5). Its source is `310_stability_direct/quasilaminar/mesh/1cyl_bodyfitted.rea`. The older generator `1cyl2d.re2` (1464 elements) is Cartesian through the disk, with nodes at r = 0 and a 0.3 pitch. Same outer domain. Regenerating from `1cyl2d.re2` does not reproduce the run mesh.

## Tests

`120_urans/check_force.py` reads the 20-time-unit moments. Pass is `|R_uu|`, `|R_vv|`, and `|R_uv|` below `1e-2` in the box `1 < x < 12`, `|y| < 2`, and centerline `u` at `x = 5` in `[0.5, 0.75]`. A missing subtraction leaves `O(1)`. On the moments now on disk the check prints `|Ruu|=6.5e-4`, `|Rvv|=1.4e-4`, `|Ruv|=1.7e-4`, `u(5,0)=0.604`, and `PASS`. The same window is `120_urans/states.png` and `120_urans/stress.png`.

```bash
python3 example/cylinder_re1m/120_urans/check_force.py
```

`RS1` is `SKIP` until `120_urans` is rebuilt and run to the last step. After that the same script fails unless `RS11cyl0.f00001` matches `E(U^2)-E(U)^2`.

The load test is a second run with `reynolds_enable` in `nekStab_usrchk` and `FRS1cyl0.f00001` present. Pass is the log line `loading frozen force FRS1cyl0.f00001` and a restart velocity that still matches the start file. A missing file must log `force not armed` and still finish. `310_stability_direct` is not this test.

The base-flow spectrum in `310_stability_direct` is a test of the imperfect SFD field only. It is not a pass for a fixed point, and it is not a mean-flow test. The mean-flow spectrum is not started.
