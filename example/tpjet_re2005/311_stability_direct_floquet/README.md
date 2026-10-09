# Triple-port jet — direct Floquet stability, Re=1900

> **Operating point: Re = 1900** (`viscosity = -1900`). See the family-level
> note in `../README.md` about the `_re2005` directory naming.

## Physics
Direct Floquet stability of the forced jet periodic orbit at Re = 1900
(St = 0.6, T = 1.66667). The forced orbit is unstable at Re = 1900 (above
Re_c ≈ 1371). The leading Floquet exponent is σ = 0.51229 with ω = ±1.23741, so the multiplier
μ = exp((σ + iω) T) is a complex pair with |μ| = 2.349 and arg μ = ±2.062 rad (118°). A pair at 118° is a secondary Hopf
(Neimark-Sacker) instability of the forced orbit. It is not a period doubling, which needs μ = −1 (arg μ = π).
Row 3 of `Spectre_NSd_conv.dat` is a real multiplier at σ = 0.0436 (μ = 1.075), also unstable.

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
sbatch run.local.slurm     # 16 ranks
```

## Result (2026-10-09, Jean Zay, `ifort`, 16 ranks)
- `Spectre_NSd_conv.dat` has ten rows. Row 1: σ = 0.5122873, ω = ±1.237406. The reference was made with 8 ranks and
  agrees to 3e-7 in σ and 1e-6 in ω.
- Wall time = 1859 s.

## Reference (ref/)
- `reference.json` — leading Floquet multiplier (σ = 0.51229, ω = 1.23741) from
  `Spectre_NSd.dat` row 1. Verify with
  `scripts/check_against_ref.py example/tpjet_re2005/311_stability_direct_floquet`.
- `ref/baseflow.png`, `ref/mode.png`, `ref/spectrum.png`.

Outputs: leading Floquet eigenmodes `dRetpjet0.f0000{1,2}` (shipped).
