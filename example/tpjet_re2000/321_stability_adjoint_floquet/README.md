# Harmonically forced jet — adjoint Floquet stability, Re = 2000

## What this stage does
The Floquet multipliers of the adjoint map of the forced orbit of `../210_baseflow_newton` (`userParam01 = 3.21`).
In exact arithmetic the adjoint map has the multipliers of the direct map of `../311_stability_direct_floquet`. Its
eigenvectors give the receptivity of the unstable multiplier to a forcing. The solver stores the orbit during the
period and runs the adjoint equations backward over it. The settings are those of the direct stage: at most
`k_dim = 128` Krylov vectors, `schur_tgt = 1` (the run stops when the leading multiplier has converged),
T = 1 / St = 1.66667. The output files are `Spectre_Ha.dat`, `Spectre_NSa_conv.dat` and the adjoint modes
`aRetpjet0.f0000N`, `aImtpjet0.f0000N`.

The mesh is axisymmetric. Nek5000 integrates the adjoint convection term on such a mesh with weights that have no
radius, and the adjoint run diverges without a correction. `nekStab_forcing` adds it (`axisymmetric_adjoint_term`;
`DOC.md`, "Adjoint on an Axisymmetric Mesh").

## Initial condition
`startFrom = BF_tpjet0.f00001`, the orbit at Re = 2000 that `../210_baseflow_newton` writes. The stage ships a
copy.

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 16 ranks
```

## Result
The first Krylov cycle converges 83 of the 128 Ritz values (tolerance 1e-6, 4492 s on 16 ranks). The leading
multiplier is real and negative, mu = -1.1573 (residual 3e-16), with growth rate sigma = ln|mu| / T = 0.0876. The
direct stage finds -1.1577, so the two differ by 4e-4 (`ref/spectrum.png`). A run with `k_dim = 40` gives -1.1579; the
scatter of 7e-4 between the two Krylov dimensions is why the reference tolerance on the leading multiplier is 2e-3.
The other converged multipliers are not those of the direct map: the second is -0.9050 against -0.8599, and -0.9654 at
`k_dim = 40`. The map is strongly non-normal, and a damped multiplier moves with the Krylov dimension, so only the
leading multiplier is compared.

The leading adjoint mode is concentrated at the lip of the nozzle, with its largest axial velocity at z = 0.64 and
r = 0.48 (`ref/mode.png`, zoomed on z < 4). The direct mode of `../311_stability_direct_floquet` has its largest
values at z = 20.6. `ref/baseflow.png` is the orbit at t = 0. `ref/reference.json` holds the numbers that
`scripts/check_against_ref.py` checks.
