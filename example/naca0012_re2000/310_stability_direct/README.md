# NACA 0012 direct global stability — Re=2000

## Physics
Direct global linear stability of the NACA 0012 steady base flow at Re = 2000
via Krylov-Schur. The leading eigenvalue is stable (σ < 0); the direct mode
localizes where a perturbation grows. Its adjoint counterpart
(`../320_stability_adjoint`) shares the same eigenvalue and localizes the
receptivity.

## nekStab Mode
- `userParam01 = 3.1` — direct stability
- `viscosity = -2000.0` — Re = 2000

## Initial condition
`startFrom = BF_naca00120.f00001` — the steady base flow from
`../210_baseflow_newton/fp`, shipped with this stage (it is the same file as in
`../320_stability_adjoint`). To regenerate it, run the 210 stage and copy its
`BF_naca00120.f00001` here.

## Run
```bash
mks naca0012
sbatch run.local.slurm     # 8 ranks
```

## Latest Result
Run of 2026-10-28 (Jean Zay, `ifort`, 8 ranks), `k_dim = 240`, `schur_tgt = 8`, `maxmodes = 10`;
wall time 2480 s (41 min 20 s, from the `logfile`).
- Leading pair: σ = -0.0665382, ω = ±8.00579 (`Spectre_NSd_conv.dat` rows 1 and 2).
- Stable pair: σ = -0.556786, ω = ±9.43325 (rows 8 and 9). An earlier run with a smaller
  target list did not list this pair; the leading pair did not change.
- Matches the adjoint eigenvalues (`../320_stability_adjoint`; σ = -0.556821, ω = ±9.43321)
  to 4e-5.

Outputs: the ten leading direct eigenmodes `dRenaca00120.f0000{1..10}` and
`dImnaca00120.f0000{1..10}` (written by the run; field files are gitignored, not tracked). The mode
file number is the row number in `Spectre_NSd_conv.dat`.

`plot.py` writes `plot_bf.png` (velocity magnitude of the base flow), `plot_spectrum.png` (converged
eigenvalues, growth rate against frequency, with the leading pair and the pair σ = -0.557 ± 9.43 i
marked) and `plot_mode.png` (vy of the leading mode in (a) and of the pair σ = -0.557 ± 9.43 i in (b)),
and `ref/mode.gif` (vy of the leading mode over one period, 24 frames). `ref/` holds copies of the
figures as `bf.png`, `spectrum.png` and `mode.png`.

## Reference (ref/)
`reference.json` holds the leading pair and the stable pair. The leading pair is read from row 1 and
the stable pair from row 9 of `Spectre_NSd_conv.dat`, because the checker cannot select a row by value.
Verify with `scripts/check_against_ref.py example/naca0012_re2000/310_stability_direct`.
