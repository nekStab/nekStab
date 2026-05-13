# Cylinder Moving Mesh — Forced Oscillation DNS

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
**Pending v2.0 validation** — to be re-run on 8 ranks and figures regenerated
through the shared `nekplot` helper during the v2.0 validation sweep. Existing `plot_*.py` scripts are legacy
matplotlib and will be consolidated into a single `plot.py` then.
