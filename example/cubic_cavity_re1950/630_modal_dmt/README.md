# cubic_cavity_re1950 / 630_modal_dmt

DMT filter bank on the limit cycle of the cubic cavity at Re = 1950 (`userParam01 = 6.4`).
The DNS runs 15 periods of the cycle with no force, and four filters watch the flow. Each filter returns the spatial structure
of the flow at one frequency as a complex field û(x), with u(x, t) ≈ Re[û(x) exp(iωt)]. Stage 610 (DMD) and stage 600 (POD) give the
same information from stored snapshots. The filter bank needs no snapshots: it runs in the time loop and writes the amplitude at the end.

## Method

A filter of angular frequency ω and width β is a damped oscillator driven by the flow u at every point:

    p' = q,    q' = u − β q − ω² p,    b = β q

`b` is u filtered around ω, with gain 1 and phase 0 at ω and gain 0 at frequency 0. The complex amplitude is
û = (b − i b'/ω) exp(−iωt), with b' from the equations, so no data are differentiated. The filter is the one of mode 1.3
(Queguineur et al. 2019), used here to analyse the flow, not to force it. `DOC.md` ("Mode 6.4") has the derivation, the leakage
of a filter into the neighbouring frequencies and the files.

## Case

- 3D flow in a cubic cavity with a uniform lid, `viscosity = -1950` (Re = 1950), 1000 elements, `lx1 = 6`, `lxd = 9`, `lpmin = 16`
  (use 16 ranks or more). No sponge, no temperature.
- `startFrom = ic_sat.f00001 time=0`: the on-orbit field of `../000_dns_period` at t = 1700, with the time reset to 0. The run
  script copies it from `../000_dns/ic_sat.f00001`, or else from `../000_dns_period/cav0.f00017`. Run `../000_dns_period` first.
- `endTime = 161.83635`: 15 periods of the cycle, T = 10.78909 (the period of `000_dns_period`). The run takes 112592 steps of 1.437e-3.
- Four filters, registered in `nekStab_usrchk` of `cav.usr` with ω₁ = 2π/T = 0.58238:

  | filter k | frequency | width β |
  |---|---|---|
  | 1 | ω₁ | 0.2 ω₁ |
  | 2 | 2 ω₁ | 0.4 ω₁ |
  | 3 | 3 ω₁ | 0.6 ω₁ |
  | 4 | 0 (the mean flow) | 0.2 ω₁ |

  The slowest filter (k = 1) settles to 1e-3 of its amplitude in 2 ln(1000)/β = 119 time units, 11 periods. The amplitude is
  averaged over the last 4 periods, from t = 11 T (`avg_from`), which also removes the ripple from the leakage of the other harmonics.
- `hist_step = 1000`: `dmt_bank.dat` has one row per 1000 steps.

## Result (40 ranks of Jean Zay, `ifort`)

Wall time 1292 s. The last field `cav0.f00001` is the flow at t = 15 T. `plot.py` rebuilds it from the filters, one harmonic at a time,

    u_end ≈ mean + Σ_{k ≤ m} Re[û_k exp(i k ω₁ t_end)],

and measures what is left, e_m = ‖u_end − mean − Σ_{k ≤ m} …‖ / ‖u_end − mean‖, in the kinetic-energy norm over the whole cavity
(Gauss-Lobatto-Legendre quadrature on the brick elements).

| k | ω | ‖û_k‖ / ‖mean‖ | amplitude at 11 T over amplitude at 15 T | e_k |
|---|---|---|---|---|
| 0 | 0 | 1 | | 1 |
| 1 | 0.5824 | 8.69e-2 | 0.9977 | 0.2274 |
| 2 | 1.1648 | 2.11e-2 | 0.9934 | 0.0591 |
| 3 | 1.7471 | 5.69e-3 | 0.9944 | 0.0192 |

- The three amplitudes change by less than 1 % over the averaging window, so the filters have settled.
- Each harmonic removes about three quarters of what is left: the fundamental, the second and the third harmonic rebuild the flow
  to 1.9 % of its fluctuation.
- The remainder is not zero. The cycle is unstable to a slow mode of period 43.5 (stage 311), and the flow has harmonics above the
  third (the DMD of stage 610 finds them up to the eighth).
- The first POD pair of stage 600 holds 94.27 % of the fluctuation energy, which leaves √(1 − 0.9427) = 0.239 of the norm after the
  fundamental. The filter bank gives 0.227 from another window of the cycle (the last 4 of 15 periods after t = 1700, against t = 1740
  to 1790 for the POD).
- The amplitude fields lie along the left wall and the bottom wall and are largest in the lower left corner, near
  (x, y) = (−0.4, −0.35), for the three harmonics; the structure gets thinner as k grows. The production of the leading Floquet
  mode (stage 411) peaks in the same corner.

`ref/dmt_bank.png`: (a) to (c) the modulus of û_k for k = 1, 2, 3 on the plane z = 0 (speed units of the lid speed U), (d) the norm of each
amplitude over its final value against t/T, with the start of the averaging window, (e) the remainder e_m against the number of harmonics.

## Reference

`ref/reference.json` checks eight quantities of a fresh run: the remainder after one, two and three harmonics (2 %, 2 % and 3 %), the norm
of each amplitude over the norm of the mean flow (1 %, 1 % and 2 %), the settling of the slowest filter (within 0.01 of 1), and the time of
the last field. A wrong frequency, filter width or period moves these numbers. `python3 ../../../scripts/check_against_ref.py .` runs the check.

## Run

```bash
mks cav
sbatch run.local.slurm        # 16 ranks or more; 1292 s on 40 ranks of Jean Zay
uv run --with numpy --with matplotlib --with pymech python plot.py
python3 ../../../scripts/check_against_ref.py .
```

## Outputs

- `bRecav0.f0000k`, `bImcav0.f0000k`: Re û_k and Im û_k of filter k, averaged over the last 4 periods. The time stamp is the end of the run.
- `dmt_bank.dat`: time, ‖u‖, ‖u − Σ b_k‖ (what no filter holds), then ‖û_k‖ for k = 1 to 4.
- `cav0.f00001`: the last field of the DNS.
- `bank_check.dat`: m, e_m, ‖û_m‖/‖mean‖. `bank_settle.dat`: k and the settling ratio. Both come from `plot.py`.
- `plot_dmt_bank.png`: copied to `ref/dmt_bank.png`.

## Literature

Queguineur et al. (2019), Phys. Fluids 31, 034101: the filter of Dynamic Mode Tracking.
