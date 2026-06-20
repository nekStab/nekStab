# Cubic Cavity Re=1950 — DNS (limit-cycle modal snapshots)

Generates the fine-sampled snapshot set consumed by POD (`600_modal_pod/`) and
DMD (`610_modal_dmd/`). Case overview: [`../README.md`](../README.md).

## Configuration
- `userParam01 = 0` (DNS), Re = 1950 (`viscosity = -1950`).
- Warm-started from the saturated-mean orbit (`ic_sat.f00001` =
  `../000_dns_period/cav0.f00066`, t = 800) so the mean flow is already settled,
  then samples ~9 periods at `writeInterval = 0.5` (~200 snapshots, double
  precision — single-precision 3D fields hang Nek MPI-IO on 16 ranks).

## Weakly-supercritical caveat (matters for DMD)
Re = 1950 sits just above the cavity's first Hopf (Re_c ≈ 1916), so the limit
cycle is very weak (centre-probe oscillation ~0.2 % of the lid speed) and
saturates **extremely slowly** — the residual transient decays at only
~7e-4 / t.u. Even by t ≈ 900 the oscillation amplitude is still drifting, so the
snapshots are not perfectly on a saturated cycle. The dominant frequency is
nonetheless sharp: the centre-probe FFT gives **St = 0.0927 (T = 10.785)**,
matching the Floquet period. The slow residual surfaces in DMD ranking — see
`../610_modal_dmd/README.md`.

## Reproduction order
`rst_cav0` → `000_dns` (initial settle) → `000_dns_period` (saturate the mean +
extract the period + Floquet orbit) → `000_dns` re-run warm from that orbit
(this stage, for better-settled modal snapshots) → `600/610`.

## Run
```bash
mks cav
sbatch run.local.slurm
```
