# Cubic Cavity Re=1950 — direct Floquet stability

Direct Floquet analysis of the 3D limit cycle. Case overview and
workflow: [`../README.md`](../README.md).

## nekStab Mode
- `userParam01 = 3.11` — direct Floquet (Krylov-Schur on the monodromy operator)
- `userParam07 = 100` — Krylov subspace dimension `k_dim`. One matvec integrates one
  period (7516 steps), about 7 min on 16 ranks, so 100 matvecs take 12.9 h.
- Period `T = 10.78909` — read from the base-flow file's time stamp
  (`BF_cav0.f00001`, stamped by `../000_dns_period/`)

## Base flow
`BF_cav0.f00001` is the on-orbit snapshot at t = 1700 of
`../000_dns_period/`, with the period T = 10.78909 written into its time
stamp by `stamp_period.py`. nekStab reads the Floquet period from that time
stamp (`matvec.f90`, `param(10) = qbase%time`). The base flow and this stage
use the same polynomial order (`lx1 = 6`); a file of another order stops the
run with a message.

## Result — an unstable pair of Floquet multipliers
69 of 100 modes converged (`eigentol = 1e-4`, residuals near 1e-16):
- **`μ = 1.0000000`** — the phase mode along the orbit. It is present
  with residual 2e-16, which confirms the period `T = 10.78909`.
- **`μ = 0.0121 ± 1.0318i`, `|μ| = 1.0319`** — a pair outside the unit circle.
  `σ = ln|μ|/T = +2.9e-3` and `ω = 0.1445`, a period of 43.5 time units, about
  four orbit periods. The limit cycle is therefore unstable to this
  oscillatory mode.
- The next multipliers are `0.9941` (real), `0.9046` (a pair) and `-0.9007`.

The adjoint run of `../321_stability_adjoint_floquet/` reproduces the pair (sigma +2.881e-3, omega 0.14457).

The pair grows slowly. The seed DNS ran 766 time units on the cycle without
visible growth, since a mode with σ = 2.9e-3 needs several thousand time units to
rise from round-off to the cycle amplitude.

Outputs: `Spectre_NSd*.dat` (σ, ω), `Spectre_Hd*.dat` (multiplier Re/Im),
`*_conv.dat` (converged-selected), and `dRe/dImcav0.f*` mode fields.
`ref/reference.json` holds the regression value of the leading pair (this
repository, not a literature value).

## Run
```bash
mks cav
sbatch run.local.slurm
python3 plot.py   # spectrum -> plot_spectrum.png
```
