# Cylinder Direct Floquet — Floquet stability of periodic orbit

## ⚠ 2D PROXY of a 3D analysis (read first)
Barkley & Henderson (1996) Mode A/B are **3D** (spanwise-periodic,
~e^{iβz}) Floquet modes of the 2D periodic wake. This case runs the
Floquet machinery in **2D only** (β = 0), so it **cannot reproduce Mode
A/B** — those *are* the third dimension; a 2D Floquet of a limit cycle
sees only the trivial unit multiplier (phase mode). This is a
**cost-reduced proxy** that exercises the nekStab pipeline at the famous
Barkley parameters, **not** a physical secondary-instability result. For
real Mode A/B, redo base flow + Floquet in 3D. Full explanation:
`../210_baseflow_newton/upo/README.md`.

## Physics
2D flow around a circular cylinder at **Re = 180** (matches the Barkley
Mode-A setup; base flow = the Re=180 UPO). Floquet analysis of the
vortex-shedding limit cycle — 2D proxy (see caveat above).

## nekStab Mode
`userParam01 = 3.11` — Direct Floquet
- `userParam07 = 100` — Krylov subspace dimension
- Sponge: left = 5, right = 5, strength = 1.7

## Prerequisites
- UPO file: `BF_1cyl0.f00001` (periodic orbit, endTime adjusted from file)

## Run
```bash
makeneks 1cyl                     # compile (links this case's .usr)
mpirun -np 8 ./nek5000 > logfile  # direct Floquet (Arnoldi on the monodromy)
```

## Results (this 2D run)
`Spectre_NSd_conv.dat` — 9 converged exponents `(σ, ω)`; multipliers
`μ = e^{(σ+iω)T}` with T = 5.1906:
- **Leading: μ = +1** (σ = −1.4×10⁻⁵ ≈ 0, ω = 0) — the trivial
  **synchronous phase mode**, exactly on the unit circle (rightmost point
  in `plot_spectrum.png`). The only neutral mode a 2D Floquet of a limit
  cycle can show.
- All other modes are **damped** (|μ| ≈ 0.94 < 1), conjugate pairs in the
  upper/lower-left of the unit disk.

NOTE: Mode A (|μ|>1 near Re~189) and Mode B (~259) are **3D** and will
**not** appear in this 2D run — see the caveat.

## Plot
```bash
# from this directory (shared renderer is example/nekplot.py)
uv run --with pymech --with scipy --with matplotlib --with numpy python plot.py
```
Outputs `plot_spectrum.png` (Floquet multipliers + unit circle),
`plot_mode.png` (leading direct mode `dRe1`, v-velocity), and `plot_bf.png`.

## Reference
Barkley & Henderson (1996), JFM 322.
