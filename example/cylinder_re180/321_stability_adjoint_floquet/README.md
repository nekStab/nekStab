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
- Sponge: left = 5, right = 5, strength = 1.7

## Prerequisites
- UPO file: `BF_1cyl0.f00001` (periodic orbit, endTime adjusted from file)

## Run
```bash
makeneks 1cyl                     # compile (links this case's .usr)
mpirun -np 8 ./nek5000 > logfile  # adjoint Floquet (Arnoldi on the adjoint monodromy)
```

## Results (this 2D run)
`Spectre_NSa_conv.dat` — converged exponents `(sigma, omega)`. The adjoint
spectrum should mirror the direct one (operator-transpose consistency). In the
2026-06-11 rerun with the sponge-strength fix, the leading adjoint Floquet
exponent is near neutral: sigma = 4.601637e-04, omega = 0. The prior saved
`run_adj.log` leading value was sigma = -1.372756e-02, omega = 0, so this is
a large shift relative to that stale output and should be treated as a real
reference update.

## Latest Result

Verified 2026-06-11 on 8 ranks (Slurm direct submission), with active sponge
strength `userParam10 = 1.7`:

- Leading adjoint Floquet exponent: sigma = 4.601637e-04, omega = 0.
  As in the direct stage, this is the trivial unit multiplier (phase
  mode, mu = 1 exactly); its deviation from zero measures the
  time-reversed orbit-replay error, larger than the direct stage's
  (-4.0e-06) but well inside the 1e-3 reference tolerance.
- Prior saved `run_adj.log` leading value: sigma = -1.372756e-02, omega = 0.
- Wall time: 1103 s.
- Outputs: `Spectre_NSa.dat`, `Spectre_Ha.dat`, `aRe1cyl0.f0000*`,
  `aIm1cyl0.f0000*`, `aRv1cyl0.f0000*`, `plot_bf.png`,
  `plot_spectrum.png`, `plot_mode.png`, and matching `ref/` copies.

## Plot
```bash
# from this directory (shared renderer is example/nekplot.py)
uv run --with pymech --with scipy --with matplotlib --with numpy python plot.py
```
Outputs `plot_spectrum.png` (adjoint Floquet multipliers + unit circle),
`plot_mode.png` (leading adjoint mode `aRe1`, v-velocity), and `plot_bf.png`.

## Reference
Barkley & Henderson (1996), JFM 322.
