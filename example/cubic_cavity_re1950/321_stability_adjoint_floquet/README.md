# Cubic Cavity Re=1950 — adjoint Floquet stability

Adjoint Floquet analysis of the 3D limit cycle, the independent check of
`../311_stability_direct_floquet/`. Case overview: [`../README.md`](../README.md).

## nekStab Mode
- `userParam01 = 3.21` — adjoint Floquet (Krylov-Schur on the adjoint monodromy operator).
  The only change against `311` is this line.
- `userParam07 = 100` — Krylov subspace dimension `k_dim`; 30 000 s (8.4 h) on 16 ranks.
- Base flow `BF_cav0.f00001`: the same stamped snapshot as in `311` (T = 10.78909).

## Result — the adjoint reproduces the direct spectrum
67 of 100 modes converged (direct: 69), residuals near 1e-16.

| multiplier | direct (311) | adjoint (this stage) |
|---|---|---|
| phase mode | 1.0000000 | 1.0000000 (`abs(mu-1) = 5e-15`) |
| unstable pair | 0.01209 ± 1.03180i, `abs(mu) = 1.03187` | 0.01133 ± 1.03151i, `abs(mu) = 1.03157` |
| sigma = ln `abs(mu)` / T | +2.908e-3 | +2.881e-3 |
| omega = arg(mu) / T | 0.14451 | 0.14457 |
| next real multiplier | 0.99414 | 0.99374 |
| next pair, `abs(mu)` | 0.90455 | 0.90592 |

The pair differs by 0.9 % in sigma and 0.05 % in omega. The cause of the
remaining difference is not analysed here. For comparison, the direct and adjoint
multipliers of `cylinder_re180` differ by 0.2 %.

Outputs: `Spectre_NSa*.dat` (sigma, omega), `Spectre_Ha*.dat` (multiplier
Re/Im), `*_conv.dat` (converged-selected), `plot_spectrum.png`.
`ref/reference.json` holds the regression value (this repository, not a
literature value).

## Run
```bash
mks cav
sbatch run.local.slurm
python3 plot.py   # spectrum -> plot_spectrum.png
```
