# Triple-port jet — direct Floquet stability, Re=1900

> **Operating point: Re = 1900** (`viscosity = -1900`). See the family-level
> note in `../README.md` about the `_re2005` directory naming.

## Physics
Direct Floquet stability of the forced jet periodic orbit at Re = 1900
(St = 0.6, T = 1.66667). The forced orbit is unstable at Re = 1900 (above
Re_c ≈ 1371); the leading Floquet multiplier μ ≈ −1 indicates a
period-doubling bifurcation.

## nekStab Mode
- `userParam01 = 3.11` — direct Floquet analysis
- `endTime = 1.666667` — orbit period (St = 0.6)
- `viscosity = -1900.0` — Re = 1900

## Initial condition
`startFrom = BF_tpjet0.f00001` — the forced orbit from
`../210_baseflow_newton` (shipped here).

## Run
```bash
mks tpjet
sbatch run.local.slurm     # 8 ranks
```

## Reference (ref/)
- `reference.json` — leading Floquet multiplier (σ = 0.51229, ω = 1.23741) from
  `Spectre_NSd.dat` row 1. Verify with
  `scripts/check_against_ref.py example/tpjet_re2005/311_stability_direct_floquet`.
- `ref/baseflow.png`, `ref/mode.png`, `ref/spectrum.png`.

Outputs: leading Floquet eigenmodes `dRetpjet0.f0000{1,2}` (shipped).
