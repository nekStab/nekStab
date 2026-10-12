# Harmonically forced jet — direct Floquet stability, Re = 2000

## What this stage does
The Floquet multipliers of the forced orbit of `../210_baseflow_newton`: the eigenvalues of the map that takes a
perturbation of the orbit to its value one forcing period later (`userParam01 = 3.11`). A multiplier outside the
unit circle is an instability of the orbit. The period is T = 1 / St = 1.66667 (St = 0.6, `endTime`). The solver
uses Krylov-Schur with at most `k_dim = 128` Krylov vectors (`userParam07`) and stops when the leading multiplier has
converged (`schur_tgt = 1` in `tpjet.usr`: the leading multiplier is real). The multipliers `mu` and the growth rates
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

## Result
The first Krylov cycle converges 81 of the 128 Ritz values (tolerance 1e-6). The orbit closes to 4.6e-5 in the
relative velocity norm. The leading multiplier is real and negative, mu = -1.1577, outside the unit circle (residual
3e-16): the forced orbit is unstable by period doubling, with growth rate sigma = ln|mu| / T = 0.0878. The next
entries of the list, three slow multipliers between 0.988 and 0.9996, had not converged after 128 steps; the
converged values below them, such as mu = -0.8599 and the pair -0.8194 +- 0.1119 i, lie inside the unit circle
(`ref/spectrum.png`). Only the leading multiplier is the stage target. The leading mode is a wave packet of the
axial velocity near the axis between z = 8 and z = 36 (`ref/mode.png`; `ref/baseflow.png` is the orbit at t = 0).
The 128 Ritz values, with the convergence flag, are in `Spectre_Hd.dat`. `ref/reference.json` holds the numbers that
`scripts/check_against_ref.py` checks. The adjoint map of `../321_stability_adjoint_floquet` has the same leading
multiplier (-1.1573 against -1.1577); its damped multipliers are not those of this map.
