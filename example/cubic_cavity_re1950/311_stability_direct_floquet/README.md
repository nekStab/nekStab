# Cubic Cavity Re=1950 — direct Floquet stability

Direct Floquet analysis of the stable 3D limit cycle. Case overview and
workflow: [`../README.md`](../README.md).

## nekStab Mode
- `userParam01 = 3.11` — direct Floquet (Krylov-Schur on the monodromy operator)
- `userParam07 = 192` — Krylov subspace dimension `k_dim`
- Period `T = 10.789` — read from the base-flow file's time stamp
  (`BF_cav0.f00001`, stamped by `../000_dns_period/`)

## Base flow
`BF_cav0.f00001` is an on-orbit DNS snapshot of the saturated limit cycle, with
its field time set to the FFT-measured period `T = 10.789`. nekStab reads the
Floquet period from that time stamp (`matvec.f90`, `param(10) = qbase%time`).

## Result — stable limit cycle
170/192 modes converged (`eigentol = 1e-3`):
- **Leading multiplier `μ = 1.000000 ± 7e-12i`** (`σ = 9.8e-13`, residual
  3.7e-13) — the trivial unit multiplier (phase/time-shift mode along the
  orbit). Its exactness confirms the period `T = 10.789`.
- **All other `|μ| < 1`**, spectrum decaying to `|μ| ≈ 3e-4` (`σ ≈ -0.75`):
  no multiplier outside the unit circle ⇒ the orbit is **stable**, consistent
  with the DNS settling onto it.

Outputs: `Spectre_NSd*.dat` (σ, ω), `Spectre_Hd*.dat` (multiplier Re/Im),
`*_conv.dat` (converged-selected), and `dRe/dImcav0.f*` mode fields.

## Run
```bash
mks cav
sbatch run.local.slurm
python3 plot.py   # spectrum -> plot_spectrum.png
```
