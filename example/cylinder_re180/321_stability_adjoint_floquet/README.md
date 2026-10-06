# Cylinder Adjoint Floquet — Adjoint Floquet stability

## ⚠ 2D PROXY of a 3D analysis (read first)
Barkley & Henderson (1996) Mode A/B are **3D** (spanwise-periodic) Floquet
modes of the 2D periodic wake. This case runs the adjoint-Floquet
machinery in **2D only** (β = 0), so it **cannot reproduce Mode A/B**
receptivity — a cost-reduced **proxy** exercising the nekStab pipeline at
the Barkley parameters, not a physical result. For real Mode A/B, redo in
3D. Full explanation: `../210_baseflow_newton/README.md`.

## Physics
2D flow around a circular cylinder at **Re = 180** (matches the Barkley
Mode-A setup; base flow = the Re=180 UPO). Adjoint Floquet analysis of
the periodic orbit for receptivity — 2D proxy (see caveat above).

## nekStab Mode
`userParam01 = 3.21` — Adjoint Floquet
- `userParam07 = 100` — Krylov subspace dimension
- Sponge off (`userParam10 = 0`). The sponge also forces the base flow, so
  with it on the orbit the Floquet run integrates is not the Newton orbit.
- Fixed `dt = 0.004` (1298 steps per period), pressure tolerance 1e-9, no filter.

## Prerequisites
- Orbit `BF_1cyl0.f00001` from `../210_baseflow_newton/`. Its time stamp is the
  period T = 5.189628.
- **`SIZE` must equal the `SIZE` of `../210_baseflow_newton/` (polynomial
  order 7).** With order 5 the orbit file is interpolated to a coarser mesh,
  the orbit no longer closes, and the multiplier of the phase mode misses 1
  by 6.8e-4 whatever the time step (measured at dt = 0.004, 0.002, 0.001).
  At order 7 the miss is 1.8e-6.

## Run
```bash
makeneks 1cyl                     # compile (links this case's .usr)
mpiexec -np 16 ./nek5000 > logfile
```
Wall time: 61 min on 16 ranks (c8).

## Results (this 2D run)
`Spectre_NSa.dat`, `Spectre_Ha_conv.dat` — multipliers `mu = exp((sigma + i omega) T)`:
- **mu = 1.0000019** (sigma = 3.8e-6, omega = 0): the trivial phase mode, as in
  the direct stage. The adjoint replays the stored orbit backwards in time.
- Four more converged multipliers: **0.8780, 0.8706, 0.8570, 0.8447** (real).
  The direct stage converged the first two as 0.8804 and 0.8724; direct and
  adjoint agree to 2.7e-3 in mu.

Mode A (|mu| > 1 near Re = 189) and Mode B (near 259) are **3D** and do not
appear in this 2D run — see the caveat above.

## Plot
```bash
# from this directory (shared renderer is example/nekplot.py)
uv run --with pymech --with scipy --with matplotlib --with numpy python plot.py
```
Outputs `plot_spectrum.png` (adjoint Floquet multipliers + unit circle),
`plot_mode.png` (leading adjoint mode `aRe1`, v-velocity), and `plot_bf.png`.

## Reference
Barkley & Henderson (1996), JFM 322.
