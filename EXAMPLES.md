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
| `example/cylinder_re100/` | 2D cylinder wake, Re = 50-100 | DNS, SFD, BoostConv, Newton, direct (linear and finite-difference) and adjoint stability, transient growth, mode animation, wavemaker, steady force sensitivity, OTD, POD/DMD/SPOD |
| `example/cylinder_re180/` | 2D cylinder wake, Re = 180 | DNS seed, Newton UPO, direct/adjoint Floquet stability, Floquet mode animation, OTD |
| `example/cylinder_re30_thermal/` | Thermally coupled cylinder wake, Re = 30 | DNS seed, Newton thermal baseflow |
| `example/cylinder_re1m/` | 2D RANS cylinder wake, Re = 1e6 | DNS, SFD base, URANS, Newton, direct stability of the SFD base (coupled and quasilaminar). Mean-flow stability not run. |
| `example/back_fstep_re500/` | Backward-facing step, Re = 500 | Newton baseflow, transient growth |
| `example/flip_flop_re62/` | Side-by-side cylinders, Re = 62 | DNS, Newton UPO, direct/adjoint Floquet stability, wavemaker |
| `example/tpjet_re2005/` | Triple-port jet, mixed Re = 1900/2005 | DMT baseflow (configuration only), forced periodic orbit, direct Floquet stability |
| `example/thermosyphon_ra500/` | Buoyancy-driven thermosyphon, Ra = 500 | DNS, Newton thermal baseflow, direct stability (linear and finite-difference) |
| `example/lid_driven_re3600/` | 2D lid-driven cavity, Re = 3600 | DNS, Newton, direct/adjoint stability |
| `example/cubic_cavity_re1914/` | 3D cubic cavity, Re = 1914 | DNS, Newton (steady), direct/adjoint stability (stable flow), OTD |
| `example/cubic_cavity_re1950/` | 3D cubic cavity, Re = 1950 | DNS, DNS period estimate, direct and adjoint Floquet stability (unstable limit cycle, pair outside the unit circle), DNS check of the unstable pair, POD/DMD |
| `example/naca0012_re2000/` | NACA 0012 airfoil, Re = 2000 | DNS, Newton, direct/adjoint stability, wavemaker |
| `example/poiseuille_re5k/` | Channel flow, Re = 5000 | OTD |
| `example/slot_fst_re495/` | Slot flow with free-stream turbulence | DNS, free-stream turbulence synthesis |
| `example/moving_cylinder_re100/` | Forced oscillating cylinder, Re = 100 | DNS forced oscillation |

## Validation status

The following examples are verified to run to completion and reproduce the
expected qualitative result with the canonical mesh and parameters shipped in
this repository.

### Verified

| Case | Re / Ra | Mode | Expected result |
|---|---:|---|---|
| `example/cylinder_re100/120_baseflow_boostconv/` | 100 | BoostConv (1.2) | converged steady state |
| `example/cylinder_re100/210_baseflow_newton/fp/` | 100 | Newton-GMRES (2.0) | converged steady state |
| `example/cylinder_re100/210_baseflow_newton/dyn/` | 100 | Newton + `ifdyntol` (2.0) | converged steady state |
| `example/cylinder_re100/110_baseflow_sfd/akervik/` | 50 | SFD baseline (1.1) | converged steady state |
| `example/cylinder_re100/110_baseflow_sfd/akervik_dyn/` | 50 | SFD + dyn-tol (1.1) | converged steady state |
| `example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs/` | 50 | SFD + OIFS (1.1) | converged steady state |
| `example/cylinder_re100/310_stability_direct/direct/` | 100 | Direct LNSE (3.1) | leading eigenvalue near sigma = 0.125 + 0.734i |
| `example/cylinder_re100/320_stability_adjoint/` | 100 | Adjoint LNSE (3.2) | matches direct spectrum |
| `example/cylinder_re100/410_postproc_animate_modes/` | 100 | Mode animation (4.50) | snapshots per period |
| `example/cylinder_re100/411_postproc_wavemaker/` | 100 | Wavemaker (4.0) | wm, sr/si, pr/pi, tr/ti fields |
| `example/cylinder_re100/412_postproc_steady_force_sensitivity/` | 100 | Steady force (4.41) | sensitivity field |
| `example/cylinder_re100/000_dns/` | 100 | DNS (0.0) | vortex shedding |
| `example/cylinder_re100/000_dns/ci_test/` | 100 | DNS smoke (0.0) | 10-step smoke |
| `example/moving_cylinder_re100/000_dns/` | 100 | DNS forced osc (0.0) | DNS subcase |
| `example/lid_driven_re3600/210_baseflow_newton/fp/` | 3600 | Newton (2.0) | converged steady state |
| `example/poiseuille_re5k/500_otd/` | 5000 | OTD smoke | builds and runs |
| `example/thermosyphon_ra500/210_baseflow_newton/` | Ra = 500 | Newton + thermal (2.0) | converged thermal steady state |
| `example/thermosyphon_ra500/310_stability_direct/direct/` | Ra = 500 | Direct + thermal (3.1) | leading real eigenvalue (pitchfork) |
| `example/naca0012_re2000/210_baseflow_newton/fp/` | 2000 | Newton (2.0) | converged steady state |
| `example/flip_flop_re62/210_baseflow_newton/` | 62 | Newton-UPO (2.1) | converged UPO |
| `example/flip_flop_re62/311_stability_direct_floquet/` | 62 | Floquet direct (3.11) | leading multiplier just above unit circle |
| `example/tpjet_re2005/210_baseflow_newton/` | 1900 | Forced UPO (2.2) | converged forced UPO |
| `example/tpjet_re2005/311_stability_direct_floquet/` | 1900 | Floquet direct (3.11) | leading multiplier above unit circle |
| `example/cubic_cavity_re1914/000_dns/` | 1914 | DNS (0.0) | steady state by t = 100 |
| `example/cubic_cavity_re1914/210_baseflow_newton/fp/` | 1914 | Newton-GMRES (2.0) | converged steady state, residual 8.4e-11 |
| `example/cubic_cavity_re1914/310_stability_direct/` | 1914 | Direct LNSE (3.1) | T1 pair sigma = -1.8e-4, omega = 0.5868 (St = 0.0934, thesis) |
| `example/cubic_cavity_re1914/320_stability_adjoint/` | 1914 | Adjoint LNSE (3.2) | matches the direct pair to 2e-4 in sigma |
| `example/cylinder_re100/330_transient_growth/` | 40 | Transient growth (3.3) | G(40) = 275; G(tau) peaks between tau = 60 and 80 |
| `example/cylinder_re180/000_dns_seed/` | 180 | DNS (0.0) | limit cycle, T = 5.1901 from the wake probe |
| `example/cylinder_re180/210_baseflow_newton/` | 180 | Newton-UPO (2.1) | converged UPO, T = 5.189628, residual 7.8e-10 (no filter) |
| `example/cylinder_re180/311_stability_direct_floquet/` | 180 | Floquet direct (3.11) | trivial multiplier 1 + 1.8e-6; next multipliers 0.880, 0.872 (real) |
| `example/cylinder_re180/321_stability_adjoint_floquet/` | 180 | Floquet adjoint (3.21) | trivial multiplier 1 + 1.9e-6; next multipliers 0.878, 0.871, 0.857, 0.845; first two match direct |

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
- `example/cubic_cavity_re1914/500_otd/` — figures come from the regularized
  lid used until 2026-10; the case file now has the uniform lid; not rerun.
- `example/tpjet_re2005/130_baseflow_dmt/` — configuration only; the DMT
  residual stalls near 0.7.
