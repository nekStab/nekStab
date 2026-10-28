# cylinder_re100 / 110_baseflow_sfd

Selective frequency damping (SFD, `userParam01 = 1.1`) computes the unstable steady wake of the circular cylinder at Re = 100 by time
marching. SFD adds the force `gain * (q - q_filtered)` to the momentum equation, where `q_filtered` is a low-pass filtered copy of
the state with the cutoff frequency `cutoff`. The damping acts on the oscillation of the unstable mode and leaves a steady state
untouched, so the time marching ends at the steady flow that Newton's method (`../210_baseflow_newton`) finds directly.

The stage has six variants in subfolders. All start from the same seed (`BF_seed_1cyl0.f00001`) and use the same mesh, 8 ranks and
`userParam04 = +-0.166`, `userParam05 = 0.26`. The sign of `userParam04` selects the formulation; the other two columns change how the
time steps are taken.

| variant | formulation | dynamic solver tolerance | time stepping | gate | converged at t | time steps | wall time on Jean Zay (s) | largest difference from Newton |
|---|---|---|---|---|---|---|---|---|
| `akervik` | Åkervik | no | BDF3, CFL 0.5 | 1e-10 | 667.9 | 66445 | 562 | 4e-05 |
| `akervik_dyn` | Åkervik | yes | BDF3, CFL 0.5 | 1e-10 | 668.1 | 66464 | 496 | 4e-05 |
| `akervik_dyn_oifs` | Åkervik | yes | OIFS, CFL 5 | 1e-10 | 765.9 | 7619 | 326 | 3e-03 |
| `casacuberta` | Casacuberta | no | BDF3, CFL 0.5 | 1e-10 | 868.9 | 86440 | 667 | 4e-05 |
| `casacuberta_dyn` | Casacuberta | yes | BDF3, CFL 0.5 | 1e-10 | 869.2 | 86462 | 588 | 4e-05 |
| `casacuberta_dyn_oifs` | Casacuberta | yes | OIFS, CFL 5 | 1e-10 | 995.3 | 9901 | 405 | 3e-03 |

Measured on Jean Zay (`ifort`, 8 ranks) on 2026-10-28. The gate is the residual (root mean square over the domain of the change of
the velocity in one time step) below which the run stops and writes `BF_1cyl0.f00001`.

## What the variants show

- The base flow does not depend on the variant. The four strict runs agree to 7 digits in every statistic of the field and differ by
  at most 4e-5 from the Newton solution (vx min -0.2001695 against -0.2001774, vx norm 264.3354 against 264.3375).
- The dynamic solver tolerance moves the time of convergence by less than 0.1 %: `akervik_dyn` stops at t = 668.13
  against 667.94, and `casacuberta_dyn` at 869.16 against 868.94. The wall time goes from 562 s to 496 s for Åkervik and
  from 667 s to 588 s for Casacuberta, about 12 % less in both formulations at the same step count (66445 against 66464,
  86440 against 86462). Wall times on Jean Zay vary between runs: in the batch of 2026-10-09 two runs of `akervik_dyn` took
  441 s and 480 s, up to 10 % apart. Read the 12 % as a trend seen in both formulations, not as a precise figure.
- With these parameters the Casacuberta filter converges more slowly than the Åkervik filter: the residual falls as
  exp(-0.0177 t) against exp(-0.0234 t) over the last fifth of the run.
- OIFS takes steps about ten times larger. At the same gate, 1e-10, it stops after 7619 steps (`akervik_dyn_oifs`, 326 s) and
  9901 steps (`casacuberta_dyn_oifs`, 405 s), against 66464 and 86462 steps for the BDF3 dyn runs. The OIFS fields stay 3.4e-3 from
  the Newton solution at any gate: at 1e-6, 1e-9 and 1e-10 the field is the same to 4e-6. The gap comes from the OIFS step at CFL 5,
  not from the gate. The strict runs stay within 4e-5.

## Parameters

`userParam04` is the frequency `f` of the oscillation to damp, in cycles per time unit, and `userParam05` its growth rate `sigma`.
Åkervik takes the cutoff `omega/2` and the gain `-2 sigma`; Casacuberta takes the cutoff `(sqrt(omega^2 + sigma^2) - sigma)/2` and the
gain `-(sqrt(omega^2 + sigma^2) + sigma)/2`, with `omega = 2 pi f`. With `f = 0.166` and `sigma = 0.26` the run prints
`cutoff 0.5215, gain -0.5200` (Åkervik) and `cutoff 0.4075, gain -0.6675` (Casacuberta).

The leading eigenvalue of the steady wake (stage `../310_stability_direct`) is 0.1248 + 0.7337 i, that is `f = 0.1168` and
`sigma = 0.1248`. The values of this stage are not those: `f = 0.166` is the shedding frequency of the saturated wake, and `sigma = 0.26`
is twice the growth rate, so the Åkervik gain is -0.52 where the rule gives -0.25. With these values every variant converges and
the residual never grows by more than 1.3 % over one time unit. For another case take `f` and `sigma` from the leading
eigenvalue of the steady flow.

## Reference and figures

Each variant has `ref/reference.json` (gate, time of convergence, four statistics of the field against the Newton values) and
`ref/*.png`. `plot_compare.png` (made by `plot_compare.py` from the residual file and `logfile` of each run) overlays the residuals of the
six variants, of BoostConv (`../120_baseflow_boostconv`, at Re = 50), of DMT (`../130_baseflow_dmt`) and of Newton
(`../210_baseflow_newton/fp`, one point per iteration) against time and against wall time. Newton writes the norm of
U(T) - q, not the change in one time step, so its curve is not the same quantity. The figure comes from the runs of
2026-10-11 (8 ranks); their wall times are not those of the table above.

Run one variant with `mks 1cyl; sbatch run.local.slurm; python3 plot.py` in its folder.

## Literature

Åkervik et al. (2006), Phys. Fluids 18, 068102; Casacuberta et al. (2018), J. Comput. Phys. 375.
