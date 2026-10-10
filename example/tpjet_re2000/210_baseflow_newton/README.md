# Harmonically forced jet — forced periodic orbit (Newton-GMRES), Re = 2000

## What this stage does
Newton-GMRES (`userParam01 = 2.2`) finds the periodic orbit of the forced jet at Re = 2000: the field q(0) for
which the flow after one forcing period, F(q), returns to q. The period is T = 1 / St = 1.66667 (St = 0.6, `endTime`).
Each Newton step solves the linear system with GMRES (`userParam07 = 128` Arnoldi vectors) on the linearized
Navier-Stokes equations along the orbit. Each linear solve stops at one fifth of the residual at the start of its
iteration (`ifdyntol` in `tpjet.usr`). The residual is the squared norm of F(q) - q (`residu_newton.dat`,
column 7) and the tolerance is 1e-8. The unstable multiplier of the orbit is not a problem for Newton's method.

## Initial condition
`startFrom = BF_Re1900_tpjet0.f00001` is the orbit of the jet at Re = 1900 (AMR_Krylov paper, polynomial order 5,
the order of this run). The Reynolds number of the run is 2000, so the first residual is not zero. The stage does
not ship the orbit that Newton writes (`BF_tpjet0.f00001`): a shipped orbit would let a run stop after one
iteration. The stages `../311_stability_direct_floquet` and `../321_stability_adjoint_floquet` ship their own copy.

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 16 ranks
```
