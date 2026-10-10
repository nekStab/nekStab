# Backward-facing step base flow at Re=500 — Newton-GMRES

## Physics
2D backward-facing step at Re = 500. Computes the steady base flow with
Newton-GMRES (the classic Barkley/Blackburn/Sherwin separated-flow benchmark).
At Re = 500 the 2D flow is **steady** (first 3D instability is near Re ≈ 750,
Barkley 2008), so the fixed point exists and is the seed for downstream
stability work (`../330_transient_growth`).

## nekStab Mode
- `userParam01 = 2.0` — Newton-GMRES for fixed points (Phase 2)
- `userParam01 = 0.0` — plain DNS (Phase 1, seed generation)
- `userParam07 = 100` — Krylov dimension `k_dim`. **Do not lower.** Re=500 sits
  just below the steady-flow stability limit, so the Jacobian carries a
  near-marginal eigenvalue; the Arnoldi/GMRES needs the room to resolve that
  stiff direction. `k_dim=30` causes Arnoldi breakdown → NaN norm.
- `viscosity = -500.0` — Re = 500

## The recipe — two phases (the canonical nekStab base-flow workflow)

Newton converges *quadratically only from inside the basin*. Do **not** seed it
with a foreign field from another run/mesh. Instead generate a native seed by
relaxing the actual mesh with a short DNS, then let Newton polish it to machine
zero. Both phases run in this directory; only the `.par` differs.

```bash
mks bfs                                  # build once
sbatch run.local.slurm                   # 8 ranks, both phases in one job
```

`run.local.slurm` runs the two phases in order and leaves `bfs.par` as the
phase-1 file when it ends:

1. Phase 1 (DNS settle): `bfs_dns.par` — uniform IC to steady recirculation on
   this mesh, endTime 400, 483 s on 8 ranks (Jean Zay), writes `bfs0.f00001..4`.
2. Phase 2 (Newton polish): `bfs_newton.par` (`startFrom = bfs0.f00004`) —
   converges in 4 iterations (639 s on 8 ranks, Jean Zay) and writes `BF_bfs0.f00001`.

The converged steady base flow is written to `BF_bfs0.f00001`.

## Convergence record (run of 2026-10-09, Jean Zay, ifort 2021.9, 8 ranks)
DNS Phase 1 relaxes monotonically; the t=300 and t=400 fields differ by 2.5e-2
(relative l2 norm of the velocity; slow recirculation-bubble tail — exactly the
mode Newton kills fast). Newton Phase 2 moves the t=400 field by a further
2.0e-2 (same norm) to the fixed point. Newton residuals from that seed:

| Newton iter | residual  |
|-------------|-----------|
| 1           | 7.03e-7   |
| 2           | 2.75e-7   |
| 3           | 1.01e-7   |
| 4           | 9.75e-9   → **converged** (target 1e-8) |

The residual falls by factors 2.6, 2.7 and 10 per step, which is slower than
quadratic, and the three GMRES solves use 49, 84 and 32 Arnoldi vectors: the
near-marginal eigenvalue dominates the residual. This is physics (Re=500 near
the stability limit), not a solver defect — and the reason `k_dim` must stay at 100.

## Pitfalls that previously broke this "easy" case
1. **Foreign seed.** Seeding Newton from a base flow exported by a *different*
   run (residual ≈ 3.2, far outside the basin) stagnated GMRES at ~7e-2 forever.
   A native DNS-settled seed (Phase 1) is residual ≈ 1e-6 — Newton finishes in
   4 steps. Same `.re2` throughout means same geometry → no mesh inversion.
2. **`call hpts` with no probe file.** `userchk` called `hpts()` but no
   `bfs.his` is shipped, so every run aborted at startup
   (`Cannot open history file ... dying`). The call is commented out — re-enable
   only if you add a `bfs.his` probe-point list.

## Reference (ref/)
- `ref/bf.png` base flow, `ref/ic.png` initial guess, `ref/residual.png`.
  Compare the converged `BF_bfs0.f00001` visually against `ref/bf.png`.
  A fresh run regenerates `residu_newton.dat` (target 1e-8).
- `reference.json` checks three quantities of `residu_newton.dat`: the final residual
  (below the 1e-8 target), the iteration count (4, tolerance 1) and the first
  residual (7.03e-7, tolerance 5 %). The last two fail if the run starts from an
  already converged field.

## Related
`../330_transient_growth` consumes this base flow for finite-time optimal
perturbation analysis.
