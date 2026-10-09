# poiseuille_re5k / 500_otd — OTD of plane Poiseuille flow at Re = 5000

## What this case shows

Three OTD modes (`lpert = 3`, `nekstab_mode = 'otd'` in the `.usr`, the same as `userParam01 = 5`) evolve with the linearized Navier-Stokes
equations about the frozen laminar profile U = 1 − y², from a white-noise start. The box is
periodic in x with length 2π and spans y = −1 to 1; Re = 1/ν = 5000 is based on the centre-line
velocity and the half-height. The mesh has 32 elements of order 7 (`lx1 = 8`), and the time
step follows a CFL number of 0.5 (about 0.05).

The linearized operator of this flow is known, so the stage checks the OTD machinery against
exact eigenvalues:

| mode | what the OTD basis converges to | exact value |
|---|---|---|
| 1 | x-independent, u = cos(πy/2), v = 0 | −ν (π/2)² = −4.934802e-4 |
| 2, 3 | the Tollmien-Schlichting wave of wavenumber α = 1 (real and imaginary part) | real part of the Orr-Sommerfeld eigenvalue −0.00175034 ± 0.26813148 i, i.e. −1.750340e-3 |

`orr_sommerfeld.py` computes the Orr-Sommerfeld eigenvalue by Chebyshev collocation. At the
critical point Re = 5772.22, α = 1.02056 it gives a growth rate of −3e-9, and at Re = 5000 it
gives the same eigenvalue to 1e-8 with 120 and with 200 points.

## Why the eigenvalues of the pair are real

The implementation removes the internal rotation of the basis in `otd_construct_linear_operator`:
it moves the lower triangle of the reduced operator `Lr` into a skew-symmetric term. The `Lr`
that is diagonalized is therefore upper triangular, its eigenvalues are its diagonal entries (the
Rayleigh quotients of the basis vectors), and `otd_eigenvalues.dat` holds real numbers. The
basis fields of modes 2 and 3 travel with the wave. They are the real and the imaginary part of
the wave, so the two entries are two real numbers whose sum is twice the real part σ of the
pair (the trace of `Lr` does not depend on the basis). The frequency of the wave is not in
`otd_eigenvalues.dat`.

## Settings

`SIZE` has `lpert = 3`. `poiseuille_OTD.usr` sets `otd_convTol = 0`, which turns the early stop
off. The stop compares the running average of the exponents between two print steps, and that
changes little long before the basis has converged (the explanation is in
`../../cylinder_re100/500_otd/README.md`). An earlier run of this stage with two modes and the
default stopped at t = 1640, where Re λ1 was −7.4e-4 against −4.9e-4 at the end.

`endTime = 5000`. The first eigenvalue approaches its exact value as exp(−2 gap t), with the
gap 1.26e-3 between mode 1 and the pair: the error falls by a factor 12 per 1000 time units.

## Result (Jean Zay, 4 ranks, `ifort`, 698 s)

| t | Re λ1 | Re λ1 + ν (π/2)² | Re λ2 | Re λ3 | running exponent of mode 1 |
|---|---|---|---|---|---|
| 500 | -1.2935e-03 | -8.00e-04 | -1.5143e-03 | -1.9566e-03 | -5.1886e-03 |
| 1000 | -1.1110e-03 | -6.18e-04 | -1.1793e-03 | -1.9665e-03 | -3.2735e-03 |
| 2000 | -6.0470e-04 | -1.11e-04 | -1.5979e-03 | -1.8228e-03 | -2.0574e-03 |
| 3000 | -5.0328e-04 | -9.80e-06 | -1.7400e-03 | -1.7564e-03 | -1.5503e-03 |
| 4000 | -4.9428e-04 | -7.99e-07 | -1.7531e-03 | -1.7693e-03 | -1.2870e-03 |
| 5000 | -4.9354e-04 | -6.47e-08 | -1.7514e-03 | -1.7527e-03 | -1.1283e-03 |

At t = 5000 the first eigenvalue is within 1.3e-4 (relative) of −ν (π/2)², and the two entries of
the pair are 0.06 % and 0.14 % from the Orr-Sommerfeld real part. A longer run (10000 time units,
4 ranks of a workstation, `gfortran`) brings the first eigenvalue to within 2e-10 of the exact value
from t = 7440. Between t = 5000 and 10000 the entries of the pair move slowly and stay within
1.0 % of the Orr-Sommerfeld value, and their mean within 0.55 % of it. The same eigenvalues come
out on Jean Zay and on the workstation to five digits from t = 500 to 4000.

The running exponents of `otd_ftle.dat` average from t = 0 and still carry the start: they are
−5.2e-3 at t = 500 and −1.1e-3 at t = 5000 for mode 1. Read `otd_eigenvalues.dat` for this
stage. An earlier reference compared the running exponents at t = 500 with −1.511e-3 and −1.772e-3.
It did not reproduce: Jean Zay and the workstation both give −5.19e-3 and −6.49e-3 there.

## Reference

`ref/reference.json` compares the last row of `otd_eigenvalues.dat` with the exact value −ν (π/2)²
(tolerance 1e-3 relative) and with the Orr-Sommerfeld real part (tolerance 0.5 % relative, for
each of the two entries of the pair). The values are not read back from an OTD run. The
tolerance of the pair belongs to the end time t = 5000: the entries move by up to 1 % later.
`python3 plot.py` writes `plot_otd_eigenvalues.png` (eigenvalues against time and the error of
the first eigenvalue) and `plot_otd_modes.png` (the three basis fields at the last output).

## Run

```bash
mks poiseuille_OTD
sbatch run.local.slurm      # 4 ranks
python3 plot.py
python3 ../../../scripts/check_against_ref.py .
```
