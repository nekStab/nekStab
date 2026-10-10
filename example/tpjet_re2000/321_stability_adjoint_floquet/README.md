# Harmonically forced jet — adjoint Floquet stability, Re = 2000

## What this stage does
The Floquet multipliers of the adjoint map of the forced orbit of `../210_baseflow_newton` (`userParam01 = 3.21`).
The adjoint map has the same multipliers as the direct map of `../311_stability_direct_floquet`. Its eigenvectors
give the receptivity of the unstable multiplier to a forcing. The solver stores the orbit during the period
and runs the adjoint equations backward over it. The settings are those of the direct stage: `k_dim = 128`,
T = 1 / St = 1.66667. The output files are `Spectre_Ha.dat`, `Spectre_NSa_conv.dat` and the adjoint modes
`aRetpjet0.f0000N`, `aImtpjet0.f0000N`.

## Initial condition
`startFrom = BF_tpjet0.f00001`, the orbit at Re = 2000 that `../210_baseflow_newton` writes. The stage ships a
copy.

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 16 ranks
```
