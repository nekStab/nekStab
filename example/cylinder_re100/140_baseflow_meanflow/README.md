# cylinder_re100 / 140_baseflow_meanflow

Mean flow of the shedding Re = 100 cylinder wake: the time average of velocity and pressure over 25 whole periods of the shedding.
The stability stage `../315_stability_direct_meanflow` uses it as the base flow.

## Why a mean flow

The stages 110 to 210 compute the steady flow, an exact solution of the Navier-Stokes equations that is unstable at Re = 100.
Its leading eigenvalue has the growth rate σ = 0.1248 and the angular frequency ω = 0.7337, which is 30 % below the frequency of
the shedding that the DNS shows (ω = 1.0417). The time average of the shedding wake is a different flow: it is not a solution of
the Navier-Stokes equations, because the Reynolds stresses of the oscillation act on it. Its recirculation bubble is much
shorter than the one of the steady flow. Barkley (2006) found that the linear operator about this mean flow has a leading eigenvalue with
σ close to zero and with ω close to the frequency of the shedding. The stage `../315_stability_direct_meanflow` repeats that analysis.

## Method

- Mode `userParam01 = 0` (DNS). `startFrom = rst_1cyl0.f00001`, the t = 100 field of `../000_dns` (copied by `run.local.slurm`).
  The DNS runs without sponge, as the steady-flow stages do.
- `userchk` in `1cyl.usr` computes `s(t)`, the integral of the velocity component v over the domain. It oscillates at the shedding
  frequency and has one upward zero crossing per period. The average opens at the first upward crossing after `tstart = 200` and
  closes at the crossing `ncycle = 25` periods later, including the step of that crossing. The routine `avg_all` of Nek5000 does
  the averaging with `param(68)` set to 1e8 in `nekStab_usrchk`, so that it never dumps and restarts the average.
- At the closing step `userchk` writes the mean through `outpost2` as `BFA1cyl0.f00001` in single precision, writes
  `meanflow_window.dat` (opening time, closing time, number of periods, period, Strouhal number; the crossing times are interpolated
  inside the step) and ends the run. `endTime = 420` is an upper bound: a run that reaches it without closing the average stops with an error.
- The probe at (5, 0) writes `1cyl.his`, which `plot.py` draws. `run.local.slurm` writes the probe list `1cyl.his` before the run,
  because the file is not tracked.

A whole number of periods leaves no trace of the oscillation in the mean. For this wake the effect is small: a mean over 150 time units
(24.87 periods, t = 250 to 400) gave the leading eigenvalue of stage 315 within 3e-5 in σ and 5e-5 in ω of the value below.

## Result (nd, `ifx`, 8 ranks)

| quantity | value |
|---|---|
| window | t = 200.632 to 351.414, 25 periods |
| period of the shedding | 6.03142 (St = 0.16580, ω = 1.04174) |
| period from the probe at (5, 0), 25 upward crossings | 6.03154 |
| range of u in the mean flow | -0.1781 to 1.2639 |
| recirculation length from the rear of the cylinder | 1.416 (steady flow: 6.154) |
| smallest u on the axis y = 0 | -0.1781 at x = 1.2 (steady flow: -0.2002 at x = 2.8) |
| wall time | 202 s |

The period from the integral of v and the period from the probe differ by 2e-5 (relative). The 2D DNS of Marquet & Carini (2020) at
Re = 100 has ω = 1.0446 for the first harmonic, 0.3 % above this run.

## Reference

`ref/reference.json` checks the number of periods, the period and the Strouhal number of `meanflow_window.dat` (1e-4), four statistics of the
mean field `BFA1cyl0.f00001` (1e-3) and the recirculation length and smallest axis velocity of `meanflow_axis.dat` (1e-2). A run that does not
close the average, or that closes it on another number of periods, fails. `python3 ../../../scripts/check_against_ref.py .` runs the check.

`plot.py` writes `meanflow_axis.dat` and `plot_meanflow.png`: (a) the probe signal with the opening and closing times, (b) u on the axis
for the mean flow and for the steady flow of `../210_baseflow_newton/fp` when that is present, (c) u of the mean flow with the line u = 0
of both flows. `ref/meanflow.png` is the figure of the run named in `ref/reference.json`.

## Run

```bash
mks 1cyl
sbatch run.local.slurm      # 8 ranks
python3 plot.py
python3 ../../../scripts/check_against_ref.py .
```

Then run `../315_stability_direct_meanflow`.

## Literature

Barkley (2006), Europhys. Lett. 75, 750. Marquet & Carini (2020), arXiv:2008.10002.
