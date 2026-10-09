# Lid-Driven Cavity — Newton-GMRES base flow

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This 2D lid-driven-cavity case is a compact Newton/baseflow exercise and a
baseline for the broader lid-driven-cavity stability discussion. The named
thesis material focuses on 3D spanwise-aspect-ratio effects, but its problem
formulation and base-flow sections define the same driven-cavity family.
Thesis source: chapter 4, `cav:sec:problem_formulation` and
`cavsubsec:BaseFlows`.

## Physics
2D lid-driven cavity at Re = 3600 with aspect ratio 1.5 (`userParam10 = 1.5`). Computes the steady base flow using Newton-GMRES.
The lid velocity is the regularized profile u = (1-(2x)^2)^2 of Shen (1991):
zero at the two lid corners, so the adjoint problem has no corner singularity.

## nekStab Mode
`userParam01 = 2` — Newton-GMRES for fixed points
- `userParam07 = 90` — Krylov subspace dimension
- `userParam10 = 1.5` — cavity aspect ratio

## Prerequisites
- Initial guess: `BF_cav0.f00001`

## Mesh
100 elements (`lelg = 100` in `SIZE`) of order 11 (`lx1 = 12`, `lxd = 18`). Compact cavity mesh.

## Run
```bash
mks cav
mpiexec -np 8 ./nek5000      # or: sbatch run.local.slurm
python3 plot.py
```

## Workflow
1. `000_dns`: DNS from the lid profile to t = 600 (steady flow), snapshots `cav0.f00001` to
   `cav0.f00006`.
2. `210_baseflow_newton/fp`: copy `cav0.f00006` of step 1 to `BF_cav0.f00001`; Newton-GMRES
   converges in two iterations and writes the base flow `BF_cav0.f00001`.
3. `310_stability_direct` and `320_stability_adjoint`: copy the `BF_cav0.f00001` of step 2.

## Resolution
The Newton and stability stages use the unfiltered time-stepper, and the DNS runs with
`filtering = explicit`, so the steady state of the DNS is a fixed point of the unfiltered map only
to the truncation error of the mesh. The table lists the first Newton residual from the t = 600
DNS state and the eigenvalues of the stability stages (4 ranks, `ifx`) at each order. The chain on Jean Zay (`ifort`) gives the same numbers to 7 digits at `lx1 = 12`:

| `lx1` | first Newton residual | direct σ₁ | adjoint σ₁ | direct pair σ ± iω | adjoint pair σ ± iω |
|---|---|---|---|---|---|
| 8 | 2.1e-6 | −0.016528 | −0.017258 | −0.04192 ± 0.78009 i | −0.04150 ± 0.78035 i |
| 10 | 1.1e-7 | −0.016764 | −0.016830 | −0.04206 ± 0.77204 i | −0.04196 ± 0.77206 i |
| 12 | 4.4e-9 | −0.016793 | −0.016796 | −0.04198 ± 0.77165 i | −0.04197 ± 0.77165 i |

At `lx1 = 6` the first residual is 3.5e-5, Newton moved to a second, unstable steady state
(growth rate +0.05, angular frequency 1.7) and the stability stages reported it; without the
filter the DNS at that order does not settle (the probe still changes at t = 600). The leading
direct and adjoint eigenvalues agree to 0.02 % at `lx1 = 12`, 0.4 % at `lx1 = 10` and 4 % at
`lx1 = 8`. The stages use `lx1 = 12`.

## Reference
Standard lid-driven cavity benchmark.
- Shen, J. (1991). Hopf bifurcation of the unsteady regularized driven cavity
  flow. *J. Comput. Phys.* 95, 228-245.
