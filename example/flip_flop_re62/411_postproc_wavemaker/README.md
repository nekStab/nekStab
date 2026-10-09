# flip_flop_re62 / 411_postproc_wavemaker

Energy budget of the two Floquet modes of the flip-flop orbit at Re = 62 (`userParam01 = 4.11`), and the overlap of the direct and adjoint modes.
The Floquet analysis (stages 311 and 321) gives the growth rate of the mode. This stage shows where the growth comes from: it integrates over one
period of the orbit the energy that the mode takes from the shear of the base flow, the energy that viscosity removes, and the energy that the
sponge removes. The growth rate that follows from these terms must equal the Floquet exponent of stage 311. That agreement is the check of the stage.

## Method

For a perturbation q of the periodic orbit U(x, t), with the energy E = 1/2 ∫ |q|² dV, the energy equation is

    dE/dt = P − D − S

- P = −∫ q_i q_j ∂U_i/∂x_j dV is the production, the work of the perturbation against the strain of the base flow.
- D = (1/Re) ∫ ∇q : ∇q̄ dV is the viscous dissipation.
- S = ∫ s(x) |q|² dV is the damping of the sponge, the force −s q that the solver adds to the perturbation equation (`spng_st = 1.7` in `2cyl.usr`).

For a Floquet mode q = exp((σ + iω) t) û(x, t) with a periodic û, the average over one period T gives

    σ = (P̄ − D̄ − S̄) / (2 Ē)

The code marches the real part (`dRe`) and the imaginary part (`dIm`) of each mode over one period on a stored copy of the orbit, and averages the
terms with the weight exp(−2σt) dt/T. The balance does not count the energy flux through the outflow boundary, so the sponge must separate the region of interest
from the outflow. The header of `src/energy_budget.f90` has the derivation and the limits.

## Case

- 2D flow around two side-by-side cylinders, `viscosity = -62` (Re = 62), 5092 elements, `lx1 = 5`, `lxd = 7`, `lpmin = 8` (use 8 ranks or more).
- `userParam01 = 4.11`, `maxmodes = 2` in `2cyl.usr`: the leading pair of Floquet modes, which has one growth rate.
- `startFrom = BF_2cyl0.f00001`, the orbit of stage 210, shipped here and in 311 and 321.
- The period is the time stamp of `BF_2cyl0.f00001`, T = 8.815172, the same period that the Floquet stages use.
  `endTime` of the `.par` file is replaced by it, and the run prints a warning when the two differ by more than 1e-6 relative.
  An earlier version of this stage carried T = 8.73356, the starting value of the Newton iteration of stage 210, which is 0.9 % short.
- One pass over the orbit takes 1224 time steps of about 7.2e-3. The run makes two passes (`dRe`, `dIm`) for each of the two modes.

## Inputs

`BF_2cyl0.f00001` and `Spectre_NSd_conv.dat` are in the folder. The modes come from the two Floquet stages; run these commands from this folder
(the modes of the committed stages are in the repository, and the eigenvalue file of 311 exists after that stage has run):

```bash
cp ../311_stability_direct_floquet/d{Re,Im}2cyl0.f0000{1,2} .
cp ../321_stability_adjoint_floquet/a{Re,Im}2cyl0.f0000{1,2} .
cp ../311_stability_direct_floquet/Spectre_NSd_conv.dat .     # only after a new run of 311
```

The budget reads the direct modes and rows 1 and 2 of `Spectre_NSd_conv.dat` (σ, ω of each mode). The adjoint modes are only used by `plot.py`.

## Result (workstation, 8 ranks, `ifx`)

Wall time 79.5 s. The Floquet exponent of stage 311 is σ = 7.89623e-3 (ω = 0.1407652).

| quantity | mode 1 | mode 2 |
|---|---|---|
| σ from the budget (line 14 of the file) | 7.9926e-3 (1.22 % above) | 7.9926e-3 |
| σ from the energy of the marched mode, ln(E(T)/E(0))/(2T) (line 16) | 7.8874e-3 (0.11 % below) | 7.8874e-3 |
| production P̄ | 3.3564e-2 | 3.3564e-2 |
| dissipation D̄ | 2.9306e-2 | 2.9306e-2 |
| sponge damping S̄ | 2.436e-4 | 2.436e-4 |
| energy Ē | 0.50234 | 0.50234 |

The terms are in the normalisation of the code, where σ = (P̄ − D̄ − S̄)/Ē. Modes 1 and 2 are the two members of one conjugate pair.

The same stage on Jean Zay (`ifort`, 8 ranks, with the modes of the 311 and 321 runs of that day, wall time 128 s) gives P̄ = 3.3563e-2, Ē = 0.50233,
σ from the budget 7.9926e-3 and σ from the energy 7.8874e-3. The two compilers agree to 5e-5 in relative terms.

- The budget closes to 1.2 %. The exponent from the energy of the marched mode involves no space derivative and is within 0.11 %, so the mode,
  the period and the stored orbit agree. The outflow is not the cause: the largest |q| on the outflow face (x = 75) is at most 5e-4 of the
  largest |q| of the mode. The remaining 1.2 % comes from the space discretisation of the production and dissipation integrals. σ is the
  difference of two terms that are eight times larger than σ Ē, so an error of 0.15 % in one of them gives 1.2 % in σ.
- Halving the time step lowers the closure error from 1.22 % to 0.99 %.
- Without the sponge term the budget gives 8.478e-3, 7.4 % above the Floquet exponent: the sponge damping shifts σ by 4.9e-4.

## Reference

`ref/reference.json` checks seven quantities of a fresh run: σ from the budget for both modes (within 2 % of the Floquet exponent), σ from the energy of
the marched mode (within 0.5 %), the sponge damping of mode 1 (2 %), and three extrema of the stored production fields (1 %; regression values).
The first two quantities fail when a term of the balance is missing: without the sponge term σ comes out 7.4 % too high. The third compares the energy
of the marched mode after one period with its energy at the start, so it tests the mode and the period with no space derivative. `python3 ../../../scripts/check_against_ref.py .` runs the check.

## Run

```bash
mks 2cyl
sbatch --ntasks=8 --mem-per-cpu=2000 run.local.slurm
# copy the modes as in Inputs, then:
uv run --with pymech --with scipy --with matplotlib --with numpy python plot.py
python3 ../../../scripts/check_against_ref.py .
```

## Outputs

- `PKE_floquet_2cyl0.f0000{1,2}`: the budget of mode 1 and mode 2, one number per line. Lines 1 to 10: integrals of the ten budget fields
  (1 to 9 production, 10 is −D̄). Line 11: P̄. Line 12: P̄ − D̄. Line 13: Ē. Line 14: σ from the budget, (line 12 − line 15)/line 13.
  Line 15: S̄. Line 16: σ from the energy of the marched mode.
- `F01*`, `F02*`: orbit-averaged production of mode 1 and mode 2. The file number is the direction of the base-flow derivative (1: x, 2: y);
  the fields `vx` and `vy` of a file belong to the x- and y-momentum equations.
- `plot_budget_field.png`: orbit-averaged production density of mode 1, the sum of the stored components. Its integral is P̄. Red: the mode takes
  energy from the base flow. Blue: it gives energy back.
- `plot_direct_mode.png`, `plot_adjoint_mode.png`: `vy` of the real part of the direct and of the adjoint mode 1 at the phase of the base-flow file.
- `plot_wavemaker.png`: product of the velocity magnitudes of the real parts of the direct and adjoint mode 1, a proxy for the wavemaker.
- The `plot_*.png` files are copied to `ref/`.

## Literature

Carini et al. (2015), JFM 778.
