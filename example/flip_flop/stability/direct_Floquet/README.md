# Flip-Flop Direct Floquet — Floquet stability of flip-flop orbit

## Physics
2D flow around two side-by-side cylinders at Re = 62. Floquet analysis of the periodic flip-flop orbit (T = 8.734) to detect secondary instabilities.

## nekStab Mode
`userParam01 = 3.11` — Direct Floquet
- `userParam07 = 48` — Krylov subspace dimension
- `endTime = 8.73356` — orbit period

## Prerequisites
- UPO file: `BF_2cyl0.f00001` (from `../../baseflow/`)

## Run
```bash
mks 2cyl        # Compile
nekbmpi 2cyl N  # Run on N MPI ranks
```

## Expected Output
Floquet multipliers of the flip-flop periodic orbit.

## Plot Notes
- Current `plot.py` shows velocity magnitude; paper plots **vx** component (vorticity proxy, asymmetry visible)
- Use PiYG diverging colormap: ±1.5 clim for UPO base flow, ±0.5 for Floquet mode
- Plot bounds: x ∈ [−3, 20], y ∈ [−3, 3]; two cylinder patches at (0, ±0.85), radius 0.5
- Neimark-Sacker bifurcation at Re_c,2 ≈ 61.17 (flip-flop mode)
- Full paper figure requires multi-Re sweep (Re = 60, 62, 63, 67) for Floquet spectrum overlay
- Reference data: Carini et al. (2014) — `carini_H67.dat`

## Reference
Carini et al. (2015), JFM 778.
