# Triple-port jet — forced periodic orbit (Newton-GMRES), Re=1900

> **Operating point: Re = 1900** (`viscosity = -1900`). See the family-level
> note in `../../README.md` about the `_re2005` directory naming.

## Physics
Axisymmetric triple-port jet at Re = 1900 (above the critical Re_c ≈ 1371).
Computes the forced periodic orbit (harmonic forcing at St = 0.6, period
T = 1/St = 1.66667) with Newton-GMRES.

## nekStab Mode
- `userParam01 = 2.2` — Newton-GMRES for forced periodic orbits
- `endTime = 1.666667` — forcing period (St = 0.6)
- `viscosity = -1900.0` — Re = 1900

## Initial condition
`startFrom = BF_tpjet0.f00001` — the converged forced orbit (shipped).
`BF_Re1900_tpjet0.f00001` is the original Re=1900 seed used to reach it.

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — converged Newton residual (`residu_newton.dat` col 6 → 1e-8).
  Verify with
  `scripts/check_against_ref.py example/tpjet_re2005/210_baseflow_newton/forced_po`.
- `ref/bf.png`, `ref/ic.png`, `ref/residual.png`.

The orbit `BF_tpjet0.f00001` feeds `../../311_stability_direct_floquet`.
