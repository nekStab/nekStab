# Cubic Cavity Re=1950 — seed DNS and period

Computes the limit cycle of the cavity at Re = 1950 with the uniform lid, and
measures its period. The period must be accurate, because the Floquet
analysis integrates over one period and its phase-mode multiplier is exactly
1 only if the integration time equals the true period.

## Workflow
1. **Seed DNS.** `ic_seed.f00001` is the steady base flow of
   `../../cubic_cavity_re1914/210_baseflow_newton/fp/` (`BF_cav0.f00001`).
   It is not a steady state at Re = 1950, so the change in Re seeds the Hopf
   mode. The run integrates to t = 3000 and is stopped by hand after the
   limit cycle is saturated (t = 1766 here). A fine `hpts` probe at the cavity
   centre writes `cav.his`. A checkpoint is written every 100 time units.
2. **Saturation.** The probe amplitude grows at about +3.7e-3 per time unit
   and is flat from t = 1000 (standard deviation 4.5e-4, in units of the lid
   speed).
3. **Period.** `python3 fft_period.py` takes the probe signal from t = 1000
   and prints the FFT peak and a zero-crossing count. A least-squares fit of
   a mean plus three harmonics with a free frequency, over t = 1400 to 1760,
   gives the sharpest value. The fit residual falls from 9e-5 (window from
   t = 1000) to 3.7e-7 (window from t = 1400), so the window must start after
   the saturation.
4. **Stamp.** nekStab reads the Floquet period from the time stamp of the
   base-flow file. `stamp_period.py` copies a checkpoint and writes T into its
   header only:

       python3 stamp_period.py cav0.f00017 10.78909 ../311_stability_direct_floquet/BF_cav0.f00001

   `cav0.f00017` (t = 1700) is also the start of `../000_dns/`.

## Result
| method | T | St |
|---|---:|---:|
| FFT peak, t = 1000–1766 | 10.7884 | 0.09269 |
| zero-crossing count, 70 cycles | 10.7984 | 0.09261 |
| harmonic fit, t = 1400–1760 | 10.78909 | 0.09269 |

The fit value T = 10.78909 is used. The zero-crossing count is 0.1 % off
because the probe carries harmonics.

## nekStab Mode
`userParam01 = 0` — plain DNS.

## Run
```bash
mks cav
cp ../../cubic_cavity_re1914/210_baseflow_newton/fp/BF_cav0.f00001 ic_seed.f00001
sbatch run.local.slurm    # 24 h limit; about 2.5 time units per minute on 16 ranks
python3 fft_period.py
```
