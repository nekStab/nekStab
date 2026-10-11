# cubic_cavity_re1950 / 411_postproc_energy_budget

Energy budget of the leading Floquet mode of the cubic cavity at Re = 1950 (`userParam01 = 4.11`).
Stage 311 gives the growth rate of the mode, σ = 2.907859e-3 (ω = 0.1445052). This stage shows where the energy of the mode
comes from: it integrates over one period of the limit cycle the energy that the mode takes from the shear of the base flow
(production) and the energy that viscosity removes (dissipation). The cavity is closed, with no sponge and no outflow, so
the balance of the two terms is exact for the continuous problem. That makes it the test of the energy route that the
open flip-flop case cannot give.

## Method

For a perturbation q of the periodic orbit U(x, t), with the energy E = 1/2 ∫ |q|² dV, the energy equation of a closed flow is

    dE/dt = P − D

- P = −∫ q_i q_j ∂U_i/∂x_j dV is the production, the work of the perturbation against the strain of the base flow.
- D = (1/Re) ∫ ∇q : ∇q̄ dV is the viscous dissipation.

For a Floquet mode q = exp((σ + iω) t) û(x, t) with a periodic û, the average over one period T gives

    σ = (P̄ − D̄) / (2 Ē)

The code marches the real part (`dRe`) and the imaginary part (`dIm`) of the mode over one period on a stored copy of the orbit, and
averages the terms with the weight exp(−2σt) dt/T. The terms are in the normalisation of the code, where σ = (P̄ − D̄)/Ē. The header of
`src/energy_budget.f90` has the derivation and the limits.

## Case

- 3D flow in a cubic cavity with a uniform lid, `viscosity = -1950` (Re = 1950), 1000 elements, `lx1 = 6`, `lxd = 9`, `lpmin = 16` (use 16 ranks or more).
- `userParam01 = 4.11`, `maxmodes = 1` in `cav.usr`: the leading Floquet mode, which is a complex pair with one growth rate.
- `startFrom = BF_cav0.f00001`, the on-orbit snapshot of `../000_dns_period` at t = 1700 with the period in its time stamp, as in stage 311.
- The period is the time stamp of `BF_cav0.f00001`, T = 10.78909. One pass over the orbit takes 7516 time steps of about 1.44e-3.
  The run builds a stored copy of the orbit once (2.5 GB per rank on 16 ranks, 39 GB in total) and then makes two passes over it (`dRe`, `dIm`).

## Inputs

The run script writes `BF_cav0.f00001` from `../000_dns_period/cav0.f00017` and copies the mode and its exponent from stage 311
(`dRecav0.f00001`, `dImcav0.f00001` and `Spectre_NSd_conv.dat`, row 1). Run `../000_dns_period` and `../311_stability_direct_floquet` first.

## Result (16 ranks, Jean Zay, `ifort`)

Wall time 711 s. A first run on 16 ranks of a 40-core node of c8 (`ifx`, 1293 s) gave the same numbers to the digits
printed below. The Floquet exponent of stage 311 is σ = 2.907859e-3.

| quantity | value |
|---|---|
| σ from the energy of the marched mode, ln(E(T)/E(0))/(2T) (line 16) | 2.90305e-3 (0.17 % below) |
| σ from the budget (line 14) | 4.6099e-3 (58.5 % above) |
| production P̄ (line 11) | 0.23284 |
| dissipation D̄ | 0.23073 |
| sponge damping S̄ (line 15) | 0 |
| energy Ē (line 13) | 0.45662 |
| P̄/Ē, D̄/Ē | 0.50992, 0.50531 |

- The exponent from the energy of the marched mode involves no space derivative. It is within 0.17 % of the Floquet exponent, so the mode, the
  period, the stored orbit and the replay of 7516 steps agree in 3D.
- The budget does not give σ here. P̄ − D̄ = 2.105e-3, and σĒ is 1.328e-3: the balance misses by 7.8e-4, which is 0.34 % of D̄.
  The exponent is the difference of two terms that are 175 times larger than σ, so a gap of 0.34 % of D̄ becomes +58.5 % of σ.
  The flip-flop stage (`../../flip_flop_re62/411_postproc_wavemaker`, 2D, open) has a gap of the same kind, 0.17 % of D̄; it becomes 1.2 % of its σ because
  its terms are only 8 times larger than σ. A gap of a few tenths of a percent is what the space discretisation of the production and dissipation
  integrals can give at `lx1 = 6`, with the production concentrated in the wall layers (figure). No resolution study was made: it needs a new orbit
  at another polynomial order.
- Use the exponent of the marched mode for σ, and the budget for where the mode takes energy. The budget gives σ when σ is not much smaller
  than P̄/Ē.

`ref/budget.png`: (a) to (c) the production density, the sum of the nine production fields, normalised by its 99.9th percentile, on the planes z = 0, y = 0 and
x = 0. Red: the mode takes energy from the base flow. Blue: it gives energy back. (d) The three growth rates.

## Reference

`ref/reference.json` checks six quantities of a fresh run: σ from the energy of the marched mode against the Floquet exponent (0.5 %),
the gap of the balance in units of D̄ (1 %), and the regression values of P̄, D̄ and Ē (0.2 %) and of the sponge damping (0, as the cavity has no sponge).
The first two compare the run with the physics; the other four fail when a term of the budget changes. `python3 ../../../scripts/check_against_ref.py .` runs the check.

## Run

```bash
mks cav
sbatch run.local.slurm        # 16 ranks, 4 GB per rank: 711 s on Jean Zay, 1293 s on a 40-core node of c8
uv run --with numpy --with pymech --with matplotlib python plot.py
python3 ../../../scripts/check_against_ref.py .
```

## Outputs

- `PKE_floquet_cav0.f00001`: the budget, one number per line. Lines 1 to 9: integrals of the nine production fields. Line 10: −D̄. Line 11: P̄.
  Line 12: P̄ − D̄. Line 13: Ē. Line 14: σ from the budget, (line 12 − line 15)/line 13. Line 15: S̄. Line 16: σ from the energy of the marched mode.
- `F01cav0.f0000{1,2,3}`: orbit-averaged production fields. The file number is the direction of the base-flow derivative (1: x, 2: y, 3: z);
  the fields `vx`, `vy` and `vz` of a file belong to the three momentum equations.
- `budget_closure.dat`: σ of the three routes, their differences, P̄, D̄, Ē and the gap of the balance in units of D̄, from `plot.py`.
- `plot_budget.png`: copied to `ref/budget.png`.

## Literature

Carini et al. (2015), JFM 778.
