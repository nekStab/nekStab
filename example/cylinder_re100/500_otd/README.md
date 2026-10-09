# cylinder_re100 / 500_otd — OTD about the steady wake (`userParam01 = 5`)

## What this case shows

An OTD basis of two orthonormal perturbations evolves with the linearized Navier-Stokes
equations about a **frozen** base flow, the steady wake of `310_stability_direct`
(`startFrom = BF_1cyl0.f00001`; the DNS is not advanced). The Gram-Schmidt step keeps the
basis orthonormal; the reduced operator `Lr` and the finite-time exponents (FTLE) follow
the subspace of fastest growth. About a steady base flow the exponents must approach the
growth rate of the leading eigenvalue pair of the 310 stage, sigma = 0.12479
(omega = 0.7337).

`lpert = 2` in `SIZE`: the leading pair is a complex pair, so two modes hold its invariant
subspace. The next eigenvalues of 310 have real parts of -0.005 to -0.03, so they are
separated from the pair by about 0.13.

OTD about a base flow that the DNS advances is `userParam01 = 5.1`; see `../../cylinder_re180/500_otd/`.

## Start

`run.local.slurm` copies `BF_1cyl0.f00001`, `dRe1cyl0.f00001` and `dIm1cyl0.f00001` from
`../310_stability_direct/direct/`. The driver reads the last two as the initial
perturbations `r01`, `r02`: the real and imaginary parts of the leading eigenmode, an
exact start on the invariant subspace. The run script first removes `r0?1cyl0.f*` files of
an earlier run, because a run writes files with these names and the driver reads the
highest number it finds.

`OTD_SEED=0 sbatch run.local.slurm` starts from white noise instead. Set `endTime = 200`
in `1cyl.par` for that start.

## No sponge

`1cyl.usr` sets `spng_st = 0`. With the sponge on, the perturbation operator differs from
the operator of the 310 stage, which has no sponge, and the exponents do not match it.

## Result (seeded start)

Growth rate of the leading pair in `310_stability_direct`: sigma = 0.124788.

| quantity | value |
|---|---|
| eigenvalues of `Lr` at t = 0 (`otd_eigenvalues.dat`, first row) | 0.124779, 0.124776 |
| FTLE at the end of the run, t = 18.2 (`otd_ftle.dat`, last row) | 0.124777, 0.124519 |

The basis starts on the invariant subspace of the pair and stays on it. The reduced
operator has real eigenvalues (the rotation at omega = 0.7337 is carried by the basis), and
their sum over the two modes is 2 sigma. The run ends by the convergence test of the
driver (see below) before `endTime = 30`.

`ref/reference.json` compares all four values with sigma of the 310 stage, not with a value
read back from an OTD run. `python3 plot.py` writes `plot_otd.png`.

### White-noise start (`OTD_SEED=0`, `endTime = 200`)

The mean of the diagonal of `Lr` over windows of 20 time units, mode 1 and mode 2:

| window | mode 1 | mode 2 | mean |
|---|---|---|---|
| 30 to 40 | 0.1532 | 0.1485 | 0.1509 |
| 40 to 60 | 0.1432 | 0.1365 | 0.1398 |
| 60 to 80 | 0.1275 | 0.1241 | 0.1258 |
| 80 to 100 | 0.1250 | 0.1244 | 0.12466 |
| 100 to 200 | 0.12478 | 0.12458 | 0.12468 |

From random fields the basis reaches the leading pair after about 100 convective time units
and then holds sigma to 1e-3 relative. The running FTLE of `otd_ftle.dat` is still 0.106 and
0.108 at t = 200, because it averages over the transient from t = 0. Use a window mean, or
`otd_FTLEPeriod`, to read the exponent of a run that starts from noise.

## Why the other OTD defaults stay

The OTD improvement pass (five candidates) changed no default. The orthonormalization uses
`lpert (lpert + 1) / 2` global reductions every `otd_gsStep = 10` steps, about one per time
step, against tens inside the Helmholtz and pressure solves of the same step. A Gram matrix
(`k_gram_matrix`) or an adaptive cadence would save a fraction of one reduction per step and
change the rounding of the results. The exponents are a running average from the start of
the run (or from the start of each window of `otd_FTLEPeriod`); a restart starts a new
average. The convergence test stops the run when the exponents change by less than
`otd_convTol = 1e-6` between two print steps, after `otd_minSteps = 200` steps. A running
average approaches its limit as `C / t`, so its change per print step is about
`C * dt_print / t^2`, smaller than the remaining error `C / t` by the factor
`dt_print / t` (about 1/20 at `t = 20`). The test can therefore stop a run before the
subspace has converged. It is left as it is, since another tolerance changes when a run
stops, not what it computes; a run that must not stop early sets a smaller `otd_convTol` or
a larger `otd_minSteps` in `nekStab_usrchk`.

## Files

- `otd_eigenvalues.dat` — time, real parts of the eigenvalues of `Lr`
- `otd_growth_rates.dat` — time, eigenvalues of the symmetric part of `Lr`
- `otd_ftle.dat` — time, FTLE of each mode
- `otd_residuals.dat` — time, change of the FTLE per print step
- `plot_otd.png` — figure from `python3 plot.py`

## Reference

Babaee and Sapsis (2016), Proc. R. Soc. A 472: 20150779.
