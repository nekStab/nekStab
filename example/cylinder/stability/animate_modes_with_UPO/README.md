# Cylinder Animate Modes with UPO — Eigenmode reconstruction on periodic orbit

## Physics
2D flow around a circular cylinder at Re = 100. Reconstructs eigenmodes superimposed on the periodic orbit over 10 periods (T = 83.20).

## nekStab Mode
`userParam01 = 4.51` — Animate modes with UPO
- `userParam07 = 100` — Krylov subspace dimension
- `endTime = 83.20` — 10x period of leading mode

## Prerequisites
- Base flow / UPO file: `BF_1cyl0.f00001`
- Converged Floquet eigenmode files

## Run
```bash
mks 1cyl        # Compile
nekbmpi 1cyl N  # Run on N MPI ranks
```

## Expected Output
Sequence of field files showing Floquet eigenmode animation on the periodic base flow.

## Reference
Internal example.
