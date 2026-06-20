# Cubic Cavity Re=1950 — period-extraction DNS

Extracts a **sharp shedding period** for the Re=1950 limit cycle so the
downstream Floquet base flow integrates over the orbit's true period. Same
recipe as `cylinder_re180/000_dns_seed_re175/`.

## Why this stage
A Floquet analysis needs the orbit period `T` to high accuracy: the trivial
unit multiplier (the phase mode) only lands exactly on `μ = 1` when the
integration window equals the true period. A coarse guess shifts the whole
spectrum. We measure `T` by FFT of a long, saturated wake-probe signal.

## Workflow
1. **Warm-start** from the saturated 3D limit cycle (`ic_orbit.f00001` =
   `../000_dns/cav0.f00179`, `t = 100`).
2. **Run a long DNS** (`endTime = 800`, ~65 periods) with a fine `hpts` probe
   at the cavity centre (`cav.his`), logging a clean single-frequency signal.
3. **FFT the probe** with `fft_period.py` (demean, Hann window, 4× zero-pad,
   parabolic sub-bin peak refine, plus cycle-counting cross-check).
4. **Stamp** an on-orbit snapshot's field time with the measured `T` and copy
   it as `BF_cav0.f00001` into `311_stability_direct_floquet/` (nekStab reads
   the Floquet period from the base-flow file's time stamp).

## Result
Measured period **T = 10.789** (`St ≈ 0.0927`). Validated downstream: the
311 Floquet trivial multiplier comes out `μ = 1.000000` to 1e-13, confirming
the period.

## nekStab Mode
`userParam01 = 0` — plain DNS.

## Run
```bash
mks cav
sbatch run.local.slurm    # writes cav.his
python3 fft_period.py     # prints T
```
