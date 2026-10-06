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
`../210_baseflow_newton/README.md`.

## Physics
2D flow around a circular cylinder at **Re = 180** (matches the Barkley
Mode-A setup; base flow = the Re=180 UPO). Floquet analysis of the
vortex-shedding limit cycle — 2D proxy (see caveat above).

## nekStab Mode
`userParam01 = 3.11` — Direct Floquet
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
Wall time: 31 min on 16 ranks (c9).

## Results (this 2D run)
`Spectre_NSd.dat`, `Spectre_Hd_conv.dat` — multipliers `mu = exp((sigma + i omega) T)`:
- **mu = 1.0000018** (sigma = 3.4e-6, omega = 0): the trivial **phase mode**
  (time shift along the orbit). Its exact value is 1; it is the only neutral
  mode a 2D Floquet of a limit cycle can show.
- Two more converged multipliers, both real: **0.8804** and **0.8724**. The
  other estimates (`Spectre_Hd.dat`) are not converged.
- A Ritz run of 100 steps converged 3 multipliers.

Mode A (|mu| > 1 near Re = 189) and Mode B (near 259) are **3D** and do not
appear in this 2D run — see the caveat above.

## Plot
```bash
# from this directory (shared renderer is example/nekplot.py)
uv run --with pymech --with scipy --with matplotlib --with numpy python plot.py
```
Outputs `plot_spectrum.png` (Floquet multipliers + unit circle),
`plot_mode.png` (leading direct mode `dRe1`, v-velocity), and `plot_bf.png`.

## Reference
Barkley & Henderson (1996), JFM 322.
