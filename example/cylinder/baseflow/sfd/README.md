# Cylinder SFD Fixed-Tolerance Testcase

This folder is the fixed-tolerance SFD baseline for the cylinder base-flow
suite.  The dynamic and OIFS variants compare against this case.

## Physics

- 2D flow around a circular cylinder with `viscosity = -50.0`
  (`Re = 50` in Nek5000's negative-viscosity Reynolds-number convention)
- steady base flow targeted with selective frequency damping
- fixed pressure and velocity solver tolerances of `1e-9`

This testcase is intentionally Re = 50.  Running the same SFD parameters at
Re = 100 produced a residual that decreased and then grew again; that was a
test-definition error for this folder, not the expected fixed-tolerance SFD
baseline.

## Active Case Definition

- `startFrom = BFRe40_1cyl0.f00001`
- shared seed hash:
  `77de81b1fa85b639a318546e36a2ef6049501c78fed284addc4a60259f48806f`
- `endTime = 400.0`
- `userParam01 = 1.1` for SFD
- `userParam04 = 0.12`
- `userParam05 = 0.05`
- `variableDt = yes`
- `timeStepper = bdf3`
- `targetCFL = 0.5`
- `lx1 = 8`, `lxd = 12`
- pressure `residualTol = 1e-9`
- velocity `residualTol = 1e-9`

The restart clock is reset to `t = 0` in `1cyl.usr` so the residual plot starts
from the testcase origin, not from the seed-file time.

## Latest Result

Latest local verification:

- date: 2026-04-28
- ranks: 8
- convergence time: `t = 168.7707`
- initial residual: `1.199387e-5` at `t = 0.005575878`
- final residual: `9.996892e-10`
- total timesteps: 30268
- elapsed wall time: 231.736 seconds
- output: converged `BF_1cyl0.f00001`

Conclusion: the corrected Re = 50 fixed-tolerance testcase converges
monotonically to the active `1e-9` gate and writes the expected `BF_` base-flow
checkpoint.

## Saved Artifacts

- `BFRe40_1cyl0.f00001`
  Initial condition shared by the SFD comparison campaign.
- `BFRe40_1cyl.nek5000`
  Collection file for visualizing the initial condition.
- `BF_1cyl0.f00001`
  Converged fixed-tolerance SFD base flow.
- `BF_1cyl.nek5000`
  Collection file for visualizing the converged base flow.
- `residu.dat`
  SFD residual history.
- `residu.png`
  Residual-decay plot.
- `plot.png`
  Combined evidence figure from `plot.py`.

## Run

```bash
mks 1cyl
mpiexec -np 8 ./nek5000
visnek BFRe40_1cyl
visnek BF_1cyl
python3 p_residu.py
python3 plot.py
```

Use 8 MPI ranks for this local testcase.  Higher rank counts on this small mesh
can trigger corrupted `elmap` metadata in Nek field files.

## Reference

Akervik et al. (2006), Phys. Fluids 18.
