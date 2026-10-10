# Harmonically forced jet — direct Floquet stability, Re = 2000

## What this stage does
The Floquet multipliers of the forced orbit of `../210_baseflow_newton`: the eigenvalues of the map that takes a
perturbation of the orbit to its value one forcing period later (`userParam01 = 3.11`). A multiplier outside the
unit circle is an instability of the orbit. The period is T = 1 / St = 1.66667 (St = 0.6, `endTime`). The solver
uses Krylov-Schur with `k_dim = 128` (`userParam07`). The multipliers `mu` and the growth rates
`sigma = ln|mu| / T` are in `Spectre_Hd.dat` (all Ritz values with a convergence flag) and `Spectre_NSd_conv.dat`
(the converged ones, as `sigma` and `omega`). A real multiplier below -1 is a period doubling.

## Initial condition
`startFrom = BF_tpjet0.f00001`, the orbit at Re = 2000 that `../210_baseflow_newton` writes. The stage ships a
copy.

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 16 ranks
```
The solver writes the leading eigenmodes `dRetpjet0.f0000N` and `dImtpjet0.f0000N`, the leading mode first.
