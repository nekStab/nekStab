# Harmonically forced jet — time-delayed feedback, Re = 2000

## What this stage does
Time-delayed feedback (TDF, `userParam01 = 1.4`) holds the flow on its unstable periodic orbit. The solver
adds the force `g (q(t) - q(t - T))` to the momentum equation, with q the velocity field and T the period of the
forcing, T = 1 / `userParam05` = 1.66667. The gain `g = -0.04432 * 2 pi / T` is fixed in the code and negative, so the
force opposes the distance from the orbit. When the flow is on the orbit, the difference is zero and the force
vanishes, so the orbit of the Navier-Stokes equations is a state without force.

The residual is the root mean square over the domain of `q(t) - q(t - T)` (`residu.dat`, column 2; column 1 is
the time, column 3 the change of the residual per step). The run ends at `endTime = 600`. When the residual falls
below `residualTol = 1e-8`, the solver writes the converged field `BF_tpjet0.f00001`. The time step follows
`targetCFL = 0.5` and is adjusted so that a period holds a whole number of steps: the delay is then the period to
machine accuracy, which the method needs.

## Initial condition
`startFrom = BF_Re1900_tpjet0.f00001` is the orbit of the jet at Re = 1900 (AMR_Krylov paper, polynomial order 5,
the order of this run). The run is at Re = 2000, so the field is close to the orbit but not on it.

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 16 ranks
```
The solver writes a snapshot every 100 time units (`tpjet0.f0000N`) and `residu.dat`.
