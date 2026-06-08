# Cylinder Adjoint Floquet — Adjoint Floquet stability

## ⚠ 2D PROXY of a 3D analysis (read first)
Barkley & Henderson (1996) Mode A/B are **3D** (spanwise-periodic) Floquet
modes of the 2D periodic wake. This case runs the adjoint-Floquet
machinery in **2D only** (β = 0), so it **cannot reproduce Mode A/B**
receptivity — a cost-reduced **proxy** exercising the nekStab pipeline at
the Barkley parameters, not a physical result. For real Mode A/B, redo in
3D. Full explanation: `../210_baseflow_newton/upo/README.md`.

## Physics
2D flow around a circular cylinder at **Re = 180** (matches the Barkley
Mode-A setup; base flow = the Re=180 UPO). Adjoint Floquet analysis of
the periodic orbit for receptivity — 2D proxy (see caveat above).

## nekStab Mode
`userParam01 = 3.21` — Adjoint Floquet
- `userParam07 = 100` — Krylov subspace dimension
- Sponge: left = 5, right = 5, strength = 1.7

## Prerequisites
- UPO file: `BF_1cyl0.f00001` (periodic orbit, endTime adjusted from file)

## Run
```bash
makeneks 1cyl                     # compile (links this case's .usr)
mpirun -np 8 ./nek5000 > logfile  # adjoint Floquet (Arnoldi on the adjoint monodromy)
```

## Results (this 2D run)
`Spectre_NSa_conv.dat` — 6 converged exponents `(σ, ω)`. The adjoint
spectrum should mirror the direct one (operator-transpose consistency),
but here the leading adjoint mode is a **decaying real mode**
(σ = −0.0137, ω = 0, |μ| ≈ 0.93) — the **neutral μ=+1 phase mode did NOT
converge into the leading adjoint set**. This is a known convergence
difficulty for the adjoint phase mode (it is weakly observable in the
adjoint operator), *not* a physics inconsistency; the oscillatory pairs
(ω ≈ 0.40, 0.46) do match their direct counterparts. To pull μ=+1 into the
adjoint set, increase the Krylov dimension / restarts. The adjoint mode
field localizes **near and upstream of the body** (receptivity), the
expected mirror of the downstream-amplifying direct mode.

## Plot
```bash
# from this directory (shared renderer is example/nekplot.py)
uv run --with pymech --with scipy --with matplotlib --with numpy python plot.py
```
Outputs `plot_spectrum.png` (adjoint Floquet multipliers + unit circle),
`plot_mode.png` (leading adjoint mode `aRe1`, v-velocity), and `plot_bf.png`.

## Reference
Barkley & Henderson (1996), JFM 322.
