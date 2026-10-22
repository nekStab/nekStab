# Cubic Cavity Re=1950 — DNS (limit-cycle modal snapshots)

Generates the snapshot set consumed by POD (`600_modal_pod/`) and DMD
(`610_modal_dmd/`): 200 snapshots, one every 0.5 time units, from t = 1700 to
1800, in double precision. Case overview: [`../README.md`](../README.md).

## Configuration
- `userParam01 = 0` (DNS), Re = 1950 (`viscosity = -1950`), uniform lid.
- Warm start `ic_sat.f00001`, a copy of `../000_dns_period/cav0.f00017`
  (t = 1700), which lies on the saturated cycle.
- Double precision: single-precision 3D snapshots hang Nek MPI-IO on 16 ranks.

## Saturation
Re = 1950 is 1.7 % above the first Hopf bifurcation (Re_c = 1916.6). The
centre-probe signal in `../000_dns_period/` grows at about +3.7e-3 per time
unit and saturates near t = 1000, at a standard deviation of 4.5e-4 (in units
of the lid speed). From t = 1000 the amplitude stays within 5 %. The snapshots
of this stage therefore lie on the cycle, with no growing transient.

## Result
The run (16 ranks, `ifx`, 3221 s) writes `cav0.f00001` to `cav0.f00200` (t = 1700.5 to 1800) and the history of the
probe at (0.25, 0.25, 0), `cav.his`. `plot.py` fits a mean and three harmonics, with a free frequency, to `vy` at the
probe and writes `cycle_check.dat` and `plot_probe.png`:

| quantity | value |
|---|---|
| period of the fit | 10.78910, against 10.78909 of `../000_dns_period` (1.1e-6 relative) |
| rms residual of the fit | 5.1e-7 |
| mean of `vy` at the probe | −0.032123 |
| half range of `vy` at the probe | 1.414e-3 |
| half range of the last period over the first | 1.0000021 |

The amplitude is the same in every period and the fit leaves a residual of 5e-7, so the 200 snapshots lie on the
cycle, with no transient.

`ref/probe.png`: (a) the probe signal and its fit, (b) the orbit in the plane of `vx` and `vy` at the probe,
(c) the half range of `vy` in each whole period.

## Reference
`ref/reference.json` checks ten quantities of a fresh run: the five numbers above, and the time and the norms of
`vx`, `vy` and `vz` of the last snapshot, `cav0.f00200`. The period (1e-5), the residual and the amplitude ratio
(1e-3) fail when the start field is off the cycle or when the Reynolds number or the mesh is wrong.
`python3 ../../../scripts/check_against_ref.py .` runs the check.

## Reproduction order
`../000_dns_period/` (seed DNS and period) → `000_dns` (this stage) →
`../600_modal_pod/` and `../610_modal_dmd/` (they read snapshots 80 to 179
through the links `cav0.f00001..100 -> ../000_dns/cav0.f00080..179`) →
`../311_stability_direct_floquet/`.

## Run
```bash
mks cav
sbatch run.local.slurm      # 36 min on 16 ranks of Jean Zay, 54 min on a 40-core node
uv run --with numpy --with scipy --with matplotlib python plot.py
```
