# NACA 0012 adjoint global stability — Re=2000

## Physics
Adjoint global linear stability of the NACA 0012 steady base flow at Re = 2000
via Krylov-Schur. The adjoint eigenvalue matches the direct one (σ, ω agree to
the convergence tolerance), as expected; the adjoint mode localizes the receptivity.

## nekStab Mode
- `userParam01 = 3.2` — adjoint stability
- `viscosity = -2000.0` — Re = 2000

## Initial condition
`startFrom = BF_naca00120.f00001` — the steady base flow from
`../210_baseflow_newton/fp` (shipped here).

## Run
```bash
mks naca0012
sbatch run.local.slurm     # 8 ranks
```

## Latest Result
Run of 2026-10-28 (Jean Zay, `ifort`, 8 ranks), `k_dim = 240`, `schur_tgt = 8`, `maxmodes = 10`,
sponge strength 1.7; wall time 1992 s (33 min 12 s, from the `logfile`).
- Leading pair: σ = -0.0665164, ω = ±8.00577 (`Spectre_NSa.dat` rows 1 and 2).
- Stable pair: σ = -0.556821, ω = ±9.43321 (rows 8 and 9), equal to the direct pair
  (`../310_stability_direct`; σ = -0.556786, ω = ±9.43325) to 4e-5.

Outputs: the ten leading adjoint eigenmodes `aRenaca00120.f0000{1..10}` and
`aImnaca00120.f0000{1..10}` (written by the run; field files are gitignored, not tracked). The mode
file number is the row number in `Spectre_NSa_conv.dat`.

`plot.py` writes `plot_bf.png` (velocity magnitude of the base flow), `plot_spectrum.png` (converged
eigenvalues, growth rate against frequency, with the leading pair and the pair σ = -0.557 ± 9.43 i
marked) and `plot_mode.png` (vy of the leading adjoint mode in (a) and of the pair σ = -0.557 ± 9.43 i
in (b); both sit at the trailing edge), and `ref/mode.gif` (vy of the leading adjoint mode over one
period, 24 frames). `ref/` holds copies of the figures as `bf.png`, `spectrum.png` and `mode.png`.

## Reference (ref/)
`reference.json` holds the leading pair and the stable pair. The leading pair is read from row 1 and
the stable pair from row 9 of `Spectre_NSa.dat`, because the checker cannot select a row by value.
Verify with `scripts/check_against_ref.py example/naca0012_re2000/320_stability_adjoint`.
