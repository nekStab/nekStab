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

## Reproduction order
`../000_dns_period/` (seed DNS and period) → `000_dns` (this stage) →
`../600_modal_pod/` and `../610_modal_dmd/` (they read snapshots 80 to 179
through the links `cav0.f00001..100 -> ../000_dns/cav0.f00080..179`) →
`../311_stability_direct_floquet/`.

## Run
```bash
mks cav
sbatch run.local.slurm      # about 36 min on 16 ranks
```
