# Examples

## Stages

Shared rationale for the numbered example stages. Geometry READMEs link here
instead of duplicating this prose.

Sources: AMR_Krylov_V5/Main sections "Theoretical framework", "Numerical
methods", "nekStab", and "Examples"; Schuh-Frantz thesis chapter 2 sections
"Theory", "Methodology", and "Validation & Verification".

### DNS

DNS stages provide the nonlinear trajectory, restart fields, and observable
signals used to seed or validate later fixed-point, periodic-orbit, stability,
and post-processing stages. They are the reference dynamics, not an eigenvalue
solver.

### Baseflow: SFD, BoostConv, Newton

Baseflow stages compute fixed points or periodic orbits that are unstable under
plain time marching. SFD and related filters damp the dominant oscillation to
approach an unstable steady state; BoostConv accelerates residual decay; the
Newton-Krylov stages solve the fixed-point or periodic-orbit residual directly
with a matrix-free time-stepper Jacobian.

### Direct Stability

Direct stability stages apply the linearized time-stepper to a base state and
use Arnoldi or Krylov-Schur iterations to approximate the leading eigenpairs.
For steady base flows this characterizes modal growth and frequency; for
periodic base flows the analogous Floquet stage studies the monodromy map over
one period.

### Adjoint Stability

Adjoint stages solve the adjoint eigenproblem associated with a direct mode.
Together, direct and adjoint modes support receptivity, structural sensitivity,
and wavemaker analyses.

### Transient Growth

Transient-growth stages target non-modal amplification. They use direct and
adjoint propagation to identify perturbations that grow over a finite horizon,
even when the modal spectrum is asymptotically stable.

### Floquet

Floquet stages start from a periodic orbit and examine perturbation growth over
one period. Multipliers inside the unit circle indicate linear decay over one
cycle; multipliers leaving the unit circle identify secondary instabilities of
the periodic orbit.

### Postproc: Wavemaker and Sensitivity

Post-processing stages combine direct and adjoint information to localize where
feedback, base-flow changes, or steady forcing most affect an eigenvalue. The
same stage family also hosts energy-budget diagnostics for interpreting the
production and dissipation mechanisms of a mode.

### Modal: POD, DMD, SPOD

Modal stages analyze saved snapshots rather than solving a new Navier-Stokes
trajectory. POD ranks coherent structures by energy, DMD by approximately
single-frequency dynamics, and SPOD by frequency-resolved coherent content.

### OTD

OTD stages evolve an orthonormal perturbation basis with the nonlinear or
time-varying base flow. They expose finite-time instability directions and
finite-time Lyapunov exponents without first freezing the dynamics into a
single steady or periodic operator.

## Case catalog

| Case | Physics | Stages in the repository |
|------|---------|--------------------------|
| `example/cylinder_re100/` | 2D cylinder wake, Re = 50-100 | DNS, SFD, BoostConv, DMT, Newton, mean flow, direct (linear and finite-difference) and adjoint stability, stability of the mean flow, transient growth, mode animation, wavemaker, OTD, POD/DMD/SPOD |
| `example/cylinder_re180/` | 2D cylinder wake, Re = 180 | DNS seed, Newton UPO, direct/adjoint Floquet stability, Floquet mode animation, OTD |
| `example/cylinder_re30_thermal/` | Thermally coupled cylinder wake, Re = 30 | DNS seed, Newton thermal baseflow |
| `example/cylinder_re1m/` | 2D RANS cylinder wake, Re = 1e6 (experimental) | DNS, SFD base, URANS, direct stability of the SFD base (coupled and quasilaminar). Not validated: no reference results, not rerun for this release, not on the gallery. Newton from the SFD field does not converge (probe documented in `110_baseflow_sfd/README.md`). Mean-flow stability not run. |
| `example/back_fstep_re500/` | Backward-facing step, Re = 500 | Newton baseflow, transient growth |
| `example/flip_flop_re62/` | Side-by-side cylinders, Re = 62 | DNS, Newton UPO, direct/adjoint Floquet stability, wavemaker |
| `example/tpjet_re2000/` | Harmonically forced jet, Re = 2000 | Time-delayed feedback and Newton for the forced periodic orbit, direct and adjoint Floquet stability |
| `example/thermosyphon_ra500/` | Buoyancy-driven thermosyphon, Ra = 500 | DNS, Newton thermal baseflow, direct stability (linear and finite-difference) |
| `example/lid_driven_re3600/` | 2D lid-driven cavity, Re = 3600 | DNS, Newton, direct/adjoint stability |
| `example/cubic_cavity_re1914/` | 3D cubic cavity, Re = 1914 | DNS, Newton (steady), direct/adjoint stability (stable flow), OTD |
| `example/cubic_cavity_re1920/` | 3D cubic cavity, Re = 1920 (above the first Hopf bifurcation) | DNS at Re = 1900, Newton (unstable steady flow), direct stability (one unstable pair) |
| `example/cubic_cavity_re1950/` | 3D cubic cavity, Re = 1950 | DNS, DNS period estimate, direct and adjoint Floquet stability (unstable limit cycle, pair outside the unit circle), DNS check of the unstable pair, POD/DMD, DMT filter bank |
| `example/naca0012_re2000/` | NACA 0012 airfoil, Re = 2000 | DNS, Newton, direct/adjoint stability, wavemaker |
| `example/poiseuille_re5k/` | Channel flow, Re = 5000 | OTD |
| `example/slot_fst_re495/` | Slot flow with free-stream turbulence | DNS, free-stream turbulence synthesis |
| `example/moving_cylinder_re100/` | Forced oscillating cylinder, Re = 100 | DNS (static and forced), Newton for the locked orbit, direct Floquet stability, DMT filter bank |

## Validation status

Each row below is a stage whose run `scripts/check_against_ref.py` compares with the numbers of its
`ref/reference.json`, and whose figures are on the validation page
([nekstab.github.io/nekStab/validation](https://nekstab.github.io/nekStab/validation/)). The table is written by
`python validation/build_gallery.py` from the summary of each reference; a stage without a complete reference is
not in it.

### Verified

<!-- verified:start (written by validation/build_gallery.py from the summary of each ref/reference.json) -->
| Stage | Re / Ra | Mode | Result |
|---|---|---|---|
| `example/cylinder_re100/000_dns/` | Re 100 | DNS (0) | DNS of the cylinder wake at Re = 100 to t = 200. The vortex street develops from a restart file. |
| `example/cylinder_re100/110_baseflow_sfd/akervik/` | Re 100 | selective frequency damping (1.1) | Steady base flow of the cylinder wake at Re = 100 by selective frequency damping, Åkervik formulation, fixed solver tolerance. |
| `example/cylinder_re100/110_baseflow_sfd/akervik_dyn/` | Re 100 | selective frequency damping (1.1) | Steady base flow of the cylinder wake at Re = 100 by selective frequency damping, Åkervik formulation, dynamic solver tolerance. |
| `example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs/` | Re 100 | selective frequency damping (1.1) | Steady base flow of the cylinder wake at Re = 100 by selective frequency damping, Åkervik formulation, dynamic solver tolerance, OIFS time stepping at CFL 5. |
| `example/cylinder_re100/110_baseflow_sfd/casacuberta/` | Re 100 | selective frequency damping (1.1) | Steady base flow of the cylinder wake at Re = 100 by selective frequency damping, Casacuberta formulation, fixed solver tolerance. |
| `example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn/` | Re 100 | selective frequency damping (1.1) | Steady base flow of the cylinder wake at Re = 100 by selective frequency damping, Casacuberta formulation, dynamic solver tolerance. |
| `example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn_oifs/` | Re 100 | selective frequency damping (1.1) | Steady base flow of the cylinder wake at Re = 100 by selective frequency damping, Casacuberta formulation, dynamic solver tolerance, OIFS time stepping at CFL 5. |
| `example/cylinder_re100/120_baseflow_boostconv/` | Re 50 | BoostConv (1.2) | Steady base flow of the cylinder wake at Re = 50 by BoostConv (Citro et al. 2017). |
| `example/cylinder_re100/130_baseflow_dmt/` | Re 100 | dynamic mode tracking (1.3) | Steady base flow of the cylinder wake at Re = 100 by dynamic mode tracking (Queguineur et al. 2019). |
| `example/cylinder_re100/140_baseflow_meanflow/` | Re 100 | DNS (0) | Mean flow of the shedding Re = 100 cylinder wake, averaged over 25 whole periods; the base flow of stage 315. |
| `example/cylinder_re100/210_baseflow_newton/fp/` | Re 100 | Newton, fixed point (2) | Steady base flow of the cylinder wake at Re = 100 by Newton-GMRES, fixed solver tolerance. |
| `example/cylinder_re100/210_baseflow_newton/dyn/` | Re 100 | Newton, fixed point (2) | Steady base flow of the cylinder wake at Re = 100 by Newton-GMRES with Eisenstat-Walker dynamic solver tolerance. |
| `example/cylinder_re100/310_stability_direct/direct/` | Re 100 | direct eigenmodes (3.1) | Leading eigenvalues of the Re = 100 cylinder wake about its steady base flow: the unstable von Kármán pair. |
| `example/cylinder_re100/315_stability_direct_meanflow/` | Re 100 | direct eigenmodes (3.1) | Leading eigenvalue of the Re = 100 cylinder wake about its mean flow: nearly marginal, at the shedding frequency. |
| `example/cylinder_re100/320_stability_adjoint/` | Re 100 | adjoint eigenmodes (3.2) | Leading adjoint eigenvalues of the Re = 100 cylinder wake about its steady base flow. |
| `example/cylinder_re100/330_transient_growth/` | Re 40 | transient growth (3.3) | Optimal energy gain of the steady cylinder wake at Re = 40 for a sweep of time horizons. |
| `example/cylinder_re100/410_postproc_animate_modes/` | Re 100 | mode animation (4.5) | Animation of the leading direct mode of the Re = 100 cylinder wake over one period, T = 8.564. |
| `example/cylinder_re100/411_postproc_wavemaker/` | Re 100 | all post-processing (4) | Wavemaker (structural sensitivity) of the Re = 100 cylinder wake from the direct and adjoint modes. |
| `example/cylinder_re100/500_otd/` | Re 100 | OTD, frozen base flow (5) | OTD basis of two modes on the frozen steady wake at Re = 100. The growth rates are compared with the leading eigenvalue of the direct stage. |
| `example/cylinder_re100/600_modal_pod/` | Re 100 |  | POD, DMD and streaming SPOD of 100 DNS snapshots of the Re = 100 limit cycle (t = 150 to 200). |
| `example/cylinder_re100/610_modal_dmd/` | Re 100 | DMD (6.2) | DMD of rank 20 of 100 DNS snapshots of the Re = 100 limit cycle (t = 150 to 200). |
| `example/cylinder_re100/620_modal_spod/` | Re 100 | SPOD (6.3) | Streaming SPOD of 100 DNS snapshots of the Re = 100 limit cycle (t = 150 to 200). |
| `example/cylinder_re180/000_dns_seed/` | Re 180 | DNS (0) | DNS of the cylinder wake at Re = 180 from a uniform flow to t = 300. It gives the seed of the periodic-orbit stage. |
| `example/cylinder_re180/210_baseflow_newton/` | Re 180 | Newton, periodic orbit (2.1) | Periodic orbit of the cylinder wake at Re = 180 (T = 5.1896) by Newton-GMRES. |
| `example/cylinder_re180/311_stability_direct_floquet/` | Re 180 | direct Floquet (3.11) | Direct Floquet analysis of the Re = 180 periodic orbit. The 2D analysis finds the trivial phase mode, multiplier 1. |
| `example/cylinder_re180/321_stability_adjoint_floquet/` | Re 180 | adjoint Floquet (3.21) | Adjoint Floquet analysis of the Re = 180 periodic orbit. The multipliers match the direct analysis. |
| `example/cylinder_re180/410_postproc_animate_modes/upo/` | Re 180 | Floquet mode animation (4.52) | Animation of the Re = 180 periodic orbit over one period, T = 5.1896. |
| `example/cylinder_re180/500_otd/` | Re 180 | OTD, base flow advanced by the DNS (5.1) | One OTD mode on the Re = 180 periodic orbit, with the base flow advanced by the DNS. Its exponent is that of the phase mode, zero. |
| `example/cylinder_re30_thermal/000_dns_seed/` | Re 25 | DNS (0) | DNS of the heated cylinder wake at Re = 25, Ri = 0.2 from a uniform flow to t = 100. |
| `example/cylinder_re30_thermal/210_baseflow_newton/` | Re 30 | Newton, fixed point (2) | Steady state of the heated cylinder wake at Re = 30, Ri = 0.2, by Newton-GMRES with the temperature field. |
| `example/moving_cylinder_re100/000_dns/` | Re 100 | DNS (0) | DNS of the static cylinder (amplitude A = 0) at Re = 100 on the mesh of the moving-cylinder case, with sponge and filter. The reference checks the frequency and the amplitude of the late probe signal. |
| `example/moving_cylinder_re100/010_dns_forced/` | Re 100 | DNS (0) | DNS of the Re = 100 cylinder in forced transverse oscillation (A = 1.08, f = 0.1643) from rest. The wake locks to the forcing: its period equals the forcing period and the field at the same phase of the forcing stops changing. |
| `example/moving_cylinder_re100/210_baseflow_newton/` | Re 100 | Newton, forced periodic orbit (2.2) | Newton-GMRES for the forced periodic orbit of the cylinder in transverse oscillation at Re = 100 (A = 1.08, f = 0.1643): the squared residual falls from 3.7e-3 to 9e-10 in 11 iterations. |
| `example/moving_cylinder_re100/311_stability_direct_floquet/` | Re 100 | direct Floquet (3.11) | Direct Floquet analysis of the forced periodic orbit of the cylinder in transverse oscillation at Re = 100: one converged multiplier, mu = 0.92988 (real, inside the unit circle), so the lock-in is stable. |
| `example/moving_cylinder_re100/630_modal_dmt/` | Re 100 | DMT filter bank (6.4) | DMT filter bank on the forced periodic orbit of the oscillating cylinder at Re = 100: the filters return the flow at the frequency of the forcing and its two harmonics, and the three harmonics rebuild the orbit to 14.5 % of its fluctuation in the near field. |
| `example/flip_flop_re62/000_dns/` | Re 62 | DNS (0) | DNS of the flip-flop configuration, two side-by-side cylinders, at Re = 62 from a uniform flow to t = 100. |
| `example/flip_flop_re62/210_baseflow_newton/` | Re 62 | Newton, periodic orbit (2.1) | Unstable periodic orbit of two side-by-side cylinders at Re = 62 (flip-flop, T = 8.8152) by Newton-GMRES. |
| `example/flip_flop_re62/311_stability_direct_floquet/` | Re 62 | direct Floquet (3.11) | Direct Floquet analysis of the flip-flop orbit at Re = 62. The leading pair lies just outside the unit circle. |
| `example/flip_flop_re62/321_stability_adjoint_floquet/` | Re 62 | adjoint Floquet (3.21) | Adjoint Floquet analysis of the flip-flop orbit at Re = 62. The leading pair matches the direct analysis. |
| `example/flip_flop_re62/411_postproc_wavemaker/` | Re 62 | Floquet energy budget (4.11) | Perturbation kinetic energy budget of the two Floquet modes of the flip-flop orbit at Re = 62. It reproduces the Floquet exponent. |
| `example/tpjet_re2000/140_baseflow_tdf/` | Re 2000 | time-delayed feedback (1.4) | Time-delayed feedback (mode 1.4) on the harmonically forced jet at Re = 2000: starting from the Re = 1900 orbit, the residual falls as t^-2.5 to 6.5e-8 at t = 600 without reaching the tolerance 1e-8. |
| `example/tpjet_re2000/210_baseflow_newton/` | Re 2000 | Newton, forced periodic orbit (2.2) | Newton-GMRES on the harmonically forced jet at Re = 2000: from the Re = 1900 orbit, 5 Newton iterations take the squared residual from 3.5e-6 to 9.7e-9 and give the forced orbit of period 1/0.6. |
| `example/tpjet_re2000/311_stability_direct_floquet/` | Re 2000 | direct Floquet (3.11) | Direct Floquet analysis of the forced jet orbit at Re = 2000: one real multiplier, -1.1577, outside the unit circle (period doubling); the other converged multipliers lie inside. |
| `example/tpjet_re2000/321_stability_adjoint_floquet/` | Re 2000 | adjoint Floquet (3.21) | Adjoint Floquet analysis of the forced jet orbit at Re = 2000: the leading multiplier is real, -1.1573, within 4e-4 of the direct multiplier; the adjoint mode sits at the nozzle lip. |
| `example/thermosyphon_ra500/000_dns/` | Ra 500 | DNS (0) | DNS of the thermosyphon at Ra = 500, Pr = 5, from rest to t = 200. |
| `example/thermosyphon_ra500/210_baseflow_newton/` | Ra 500 | Newton, fixed point (2) | Steady base flow of the thermosyphon at Ra = 500 by Newton-GMRES. |
| `example/thermosyphon_ra500/310_stability_direct/direct/` | Ra 500 | direct eigenmodes (3.1) | Direct stability of the thermosyphon at Ra = 500. The leading eigenvalue is real and slightly positive (pitchfork). |
| `example/back_fstep_re500/210_baseflow_newton/` | Re 500 |  | Steady base flow of the backward-facing step at Re = 500: DNS to settle, then Newton-GMRES. |
| `example/back_fstep_re500/330_transient_growth/` | Re 500 | transient growth (3.3) | Optimal transient growth of the backward-facing step at Re = 500 at the horizon tau = 29.913 of Blackburn, Barkley and Sherwin (2008). |
| `example/lid_driven_re3600/000_dns/` | Re 3600 | DNS (0) | DNS of the lid-driven cavity at Re = 3600 (aspect ratio 1.5, regularized lid) to t = 600. The flow relaxes to its steady state. |
| `example/lid_driven_re3600/210_baseflow_newton/fp/` | Re 3600 | Newton, fixed point (2) | Steady base flow of the lid-driven cavity at Re = 3600 by Newton-GMRES from the DNS state. |
| `example/lid_driven_re3600/310_stability_direct/` | Re 3600 | direct eigenmodes (3.1) | Leading eigenvalues of the lid-driven cavity at Re = 3600: a stable real eigenvalue and a stable complex pair. |
| `example/lid_driven_re3600/320_stability_adjoint/` | Re 3600 | adjoint eigenmodes (3.2) | Adjoint eigenvalues of the lid-driven cavity at Re = 3600. They match the direct eigenvalues. |
| `example/cubic_cavity_re1914/000_dns/` | Re 1914 | DNS (0) | DNS of the cubic cavity at Re = 1914 from rest to t = 100. It gives the start field of the Newton stage. |
| `example/cubic_cavity_re1914/210_baseflow_newton/fp/` | Re 1914 | Newton, fixed point (2) | Steady base flow of the cubic cavity at Re = 1914 by Newton-GMRES from the DNS field. |
| `example/cubic_cavity_re1914/310_stability_direct/` | Re 1914 | direct eigenmodes (3.1) | Leading eigenpairs of the steady cubic cavity flow at Re = 1914. The leading pair is slightly stable: sigma = -1.83e-4, omega = 0.5868. |
| `example/cubic_cavity_re1914/320_stability_adjoint/` | Re 1914 | adjoint eigenmodes (3.2) | Leading adjoint eigenpair of the steady cubic cavity flow at Re = 1914. It equals the direct pair. |
| `example/cubic_cavity_re1914/330_transient_growth/` | Re 1914 | transient growth (3.3) | Optimal transient growth of the steady cubic-cavity flow at Re = 1914 for the horizon tau = 10, with the optimal perturbation and its response. |
| `example/cubic_cavity_re1914/411_postproc_wavemaker/` | Re 1914 | all post-processing (4) | Wavemaker (structural sensitivity) of the cubic cavity at Re = 1914 from the direct and adjoint modes. |
| `example/cubic_cavity_re1914/500_otd/` | Re 1914 | OTD, frozen base flow (5) | Two OTD modes on the steady flow of the cubic lid-driven cavity at Re = 1914, started on the eigenplane of the leading pair of 310_stability_direct. The basis must stay on that plane, and the mean of the two diagonal entries of the reduced operator must equal the growth rate of the pair. |
| `example/cubic_cavity_re1920/000_dns/` | Re 1900 | DNS (0) | DNS of the cubic cavity at Re = 1900 from rest to t = 100. It gives the start field of the Newton stage at Re = 1920. |
| `example/cubic_cavity_re1920/210_baseflow_newton/fp/` | Re 1920 | Newton, fixed point (2) | Steady base flow of the cubic cavity at Re = 1920, above the critical Reynolds number, by Newton-GMRES. |
| `example/cubic_cavity_re1920/310_stability_direct/` | Re 1920 | direct eigenmodes (3.1) | Leading eigenpairs of the steady cubic cavity flow at Re = 1920. One pair is unstable: sigma = +5.84e-4, omega = 0.5868. |
| `example/cubic_cavity_re1950/000_dns/` | Re 1950 | DNS (0) | DNS of the cubic cavity at Re = 1950 on the saturated limit cycle from t = 1700 to 1800. It writes the 200 snapshots that the POD and DMD stages read. |
| `example/cubic_cavity_re1950/000_dns_period/` | Re 1950 | DNS (0) | DNS of the cubic cavity at Re = 1950 until the limit cycle has saturated, and the period of the cycle from the probe signal. |
| `example/cubic_cavity_re1950/311_stability_direct_floquet/` | Re 1950 | direct Floquet (3.11) | Direct Floquet analysis of the limit cycle of the cubic cavity at Re = 1950. One pair lies outside the unit circle, \|mu\| = 1.032. |
| `example/cubic_cavity_re1950/321_stability_adjoint_floquet/` | Re 1950 | adjoint Floquet (3.21) | Adjoint Floquet analysis of the limit cycle of the cubic cavity at Re = 1950. The pair matches the direct analysis. |
| `example/cubic_cavity_re1950/350_dns_unstable_orbit/` | Re 1950 | DNS (0) | DNS of the unstable limit cycle of the cubic cavity at Re = 1950 with a small part of the leading Floquet mode. The growth reproduces the Floquet pair. |
| `example/cubic_cavity_re1950/411_postproc_energy_budget/` | Re 1950 | Floquet energy budget (4.11) | Energy budget of the leading Floquet mode of the cubic cavity at Re = 1950. The energy of the marched mode gives the Floquet exponent to 0.2 %; production and dissipation balance to 0.34 % of the dissipation. |
| `example/cubic_cavity_re1950/600_modal_pod/` | Re 1950 | POD, DMD and SPOD (6) | POD of 100 snapshots of the limit cycle of the cubic cavity at Re = 1950. Two modes hold 94.3 % of the energy and turn at the frequency of the cycle. |
| `example/cubic_cavity_re1950/610_modal_dmd/` | Re 1950 | POD, DMD and SPOD (6) | DMD of 100 snapshots of the limit cycle of the cubic cavity at Re = 1950. The modes lie on the harmonics of the cycle, 1 to 6 times 2 pi / T, and none grows. |
| `example/cubic_cavity_re1950/630_modal_dmt/` | Re 1950 | DMT filter bank (6.4) | DMT filter bank on the limit cycle of the cubic cavity at Re = 1950: the filters return the spatial structure of the flow at the frequency of the cycle and its two harmonics, and the three harmonics rebuild the flow to 1.9 % of its fluctuation. |
| `example/naca0012_re2000/000_dns/` | Re 2000 | DNS (0) | DNS of the NACA 0012 wake at Re = 2000. The wake rings down at the frequency and rate of the leading eigenvalue. |
| `example/naca0012_re2000/210_baseflow_newton/fp/` | Re 2000 | Newton, fixed point (2) | Steady base flow of the NACA 0012 airfoil at Re = 2000 by Newton-GMRES from a DNS snapshot. |
| `example/naca0012_re2000/310_stability_direct/` | Re 2000 | direct eigenmodes (3.1) | Leading eigenvalues of the NACA 0012 steady flow at Re = 2000. The flow is stable. |
| `example/naca0012_re2000/320_stability_adjoint/` | Re 2000 | adjoint eigenmodes (3.2) | Leading adjoint eigenvalues of the NACA 0012 steady flow at Re = 2000. |
| `example/naca0012_re2000/411_postproc_wavemaker/` | Re 2000 | all post-processing (4) | Wavemaker, base-flow sensitivity and energy budget of the NACA 0012 modes at Re = 2000. |
| `example/poiseuille_re5k/500_otd/` | Re 5000 |  | Three OTD modes on the laminar channel flow at Re = 5000. The growth rates are compared with theory: viscous decay and the Orr-Sommerfeld eigenvalue. |
| `example/slot_fst_re495/000_dns/` | Re 495 | DNS (0) | DNS of a slot jet at Re = 495 with a synthetic free-stream turbulence inflow, to t = 100. |
<!-- verified:end -->

### Reproducing

Each validated case includes a `README.md` describing physics, mode selection,
and expected output, plus the necessary `.par`, `.usr`, `SIZE`, and where
relevant, base-flow files.

```bash
cd example/<path>
mks <casename>
nekbmpi <casename> <N>     # or sbatch <case>.slurm on a cluster
python3 plot.py            # regenerate the evidence figure
```

### Not yet verified

Stages not in the table above ship their inputs, and their `README.md` gives
their state. Known gaps:

- `example/cylinder_re30_thermal/` — baseflow only, no stability stage yet.
- `example/cylinder_re1m/` — experimental, see the case catalog.
