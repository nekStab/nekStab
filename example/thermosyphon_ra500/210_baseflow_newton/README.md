# Thermosyphon base flow at Ra=500 — Newton-GMRES

## Physics
Closed-loop thermosyphon (natural-convection thermal cavity) at Rayleigh
number Ra = 500, slightly above the critical Ra_c ≈ 494. Computes the steady
base flow with Newton-GMRES.

## nekStab Mode
- `userParam01 = 2.0` — Newton-GMRES for fixed points
- `userParam06 = 500.0` — Rayleigh number Ra
- Buoyancy-driven (Boussinesq) thermal case.

## Initial condition
`startFrom = BF_Ra400_tsyphon0.f00001` — a nearby Ra=400 base flow as the
Newton initial guess (shipped). Converges to `BF_tsyphon0.f00001` at Ra=500.

## Run
```bash
mks tsyphon
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — converged Newton residual (`residu_newton.dat` col 6, last
  iteration → ~1.5e-11). Verify with
  `scripts/check_against_ref.py example/thermosyphon_ra500/210_baseflow_newton`.
- `ref/bf.png` base flow, `ref/ic.png` initial guess, `ref/residual.png` decay.

The converged `BF_tsyphon0.f00001` feeds `../310_stability_direct/direct`.
