# Cylinder Moving Mesh — Forced Oscillation DNS

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This case exercises forced motion and moving-domain mechanics for the circular
cylinder wake. The named thesis chapters discuss periodic orbits, Floquet
multipliers, and circular-cylinder wake instabilities, but they do not contain
a separate forced-oscillating-cylinder case section. Thesis source for the
underlying concepts: chapter 2, `Periodic orbits`, `Bifurcations of periodic
orbits`, and `sec:ex:1cyl`.

## Physics
2D flow around a circular cylinder undergoing forced transverse oscillation
at amplitude `A = 1.08` (`userParam05`) and frequency `f = 0.1643` (`userParam06`),
driven by the `my_meshv` user routine in `1cyl.usr`. DNS at the underlying
flow regime; demonstrates the moving-mesh / forced-oscillation capability
of nekStab.

## nekStab Mode
`userParam01 = 0` — standard DNS

Alternate modes available in the same case file:
- `userParam01 = 2.2` — Newton for forced UPO (uncomment in `.par`)
- `userParam01 = 3.11` — Floquet direct eigenproblem (uncomment in `.par`)

`userParam07 = 200` — Krylov subspace dimension (used when running Newton/Floquet variants).

## Prerequisites
- Base flow: `BF_1cyl0.f00001` (saved seed, included)

## Run
```bash
mks 1cyl
nekbmpi 1cyl 8
```

## Expected Output
DNS time-history with the cylinder oscillating about its mean position; the
flow phase-locks to the forcing once the transient decays.

## Status
**DEFERRED** (2026-05-21) — the whole `moving_cylinder_re100/` family is
deferred until the baseline Re=100 ladder (`cylinder_re100/` and
`cylinder_re180/` Floquet) is fully validated. The 4 catalog cases
under `flow_family='moving_cylinder'` (`newton_forced_upo`,
`floquet_direct`, `energy_response`, plus the legacy DNS entry) all
remain at `status=DEFERRED`. No active work; revisit when the static
ladder is complete.

When picking up the work, plan to:

- re-run the DNS at amp=1.08 / f=0.1643 on 8 ranks
- regenerate figures through the shared `nekplot` helper
- consolidate the legacy `plot_*.py` scripts into a single `plot.py`
- compute Newton-forced-UPO (uparam=2.2) using the prescribed period
  `T = 1/0.1643 ≈ 6.09` as the initial guess
- once the UPO is converged, run Floquet direct (uparam=3.11) to
  characterize the period-T orbit's 3D-mode response
