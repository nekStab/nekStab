# Cubic Cavity Re=1950 — seed DNS and period

Computes the limit cycle of the cavity at Re = 1950 with the uniform lid, and
measures its period. The period must be accurate, because the Floquet
analysis integrates over one period and its phase-mode multiplier is exactly
1 only if the integration time equals the true period.

## Workflow
1. **Seed DNS.** `ic_seed.f00001` is the steady base flow of
   `../../cubic_cavity_re1914/210_baseflow_newton/fp/` (`BF_cav0.f00001`).
   It is not a steady state at Re = 1950, so the change in Re seeds the Hopf
   mode. `cav.par` sets `endTime = 3000.0`, but the period needs the run only
   to t = 1800: set `endTime = 1800.0`, or stop the job after `cav0.f00018` is
   written. The validation run used `endTime = 1800.0` (11 h on 16 ranks). A
   fine `hpts` probe at the cavity centre writes `cav.his`. A checkpoint is
   written every 100 time units.
2. **Saturation.** The oscillation amplitude of the probe grows at +3.7e-3 per
   time unit (panel (a) of `plot_period.png`). It overshoots its final value by
   2 % near t = 1100 and is within 0.05 % of it from t = 1300. The standard
   deviation of the saturated signal is 4.5e-4, in units of the lid speed.
3. **Period.** `python3 plot.py` reads `cav.his` and measures the period over
   t = 1400 to the end of the run with three methods: a least-squares fit of a
   mean plus three harmonics with a free frequency (the value used), the peak
   of the FFT, and the mean time between zero crossings. It writes
   `period_fit.dat` and `plot_period.png`. The rms residual of the fit, in
   units of the lid speed, is 3.6e-7 for a window from t = 1400, 3.0e-6 from
   t = 1200 and 9.2e-5 from t = 1000, so the window must start after the
   saturation. The older `fft_period.py` prints the FFT peak and the crossing
   count from t = 1000.
4. **Stamp.** nekStab reads the Floquet period from the time stamp of the
   base-flow file. `stamp_period.py` copies a checkpoint and writes T into its
   header only:

       python3 stamp_period.py cav0.f00017 10.78909 ../311_stability_direct_floquet/BF_cav0.f00001

   `cav0.f00017` (t = 1700) is also the start of `../000_dns/`.

## Result
Window t = 1400–1800 (`period_fit.dat`):

| method | T | St |
|---|---:|---:|
| harmonic fit | 10.78909 | 0.09269 |
| FFT peak | 10.78915 | 0.09269 |
| zero-crossing count | 10.78910 | 0.09269 |

The fit value T = 10.78909 is used. The three methods agree to 1e-5. The same
methods on a window from t = 1000 give 10.7870 (fit), 10.7885 (FFT) and 10.7980
(zero crossings): that window still holds the end of the saturation. The probe
at the cavity centre has the mean 0.026736 and the half-range 6.67e-4, in units
of the lid speed.

## nekStab Mode
`userParam01 = 0` — plain DNS.

## Run
```bash
mks cav
cp ../../cubic_cavity_re1914/210_baseflow_newton/fp/BF_cav0.f00001 ic_seed.f00001
sbatch run.local.slurm    # 24 h limit; about 2.7 time units per minute on 16 ranks
python3 plot.py
```

## Reference
`ref/reference.json` checks the three estimates of the period (1e-5 for the fit, 1e-4 for the FFT peak,
2e-5 for the zero crossings), the rms residual of the fit (below 1.4e-6), the mean of the probe signal
(1e-3) and its half-range (2e-3). `plot.py` writes `period_fit.dat`, so the check fails on a clean
checkout until the run has been made.

## Figure
`python3 plot.py` writes `plot_period.png`: (a) the maximum and the minimum of the probe signal in each
period, with the start of the fit window, and (b) the last five periods of the signal against the
harmonic fit. The copy `ref/period.png` is the figure of the run named in `ref/reference.json`.
