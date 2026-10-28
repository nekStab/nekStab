# cylinder_re100 / 315_stability_direct_meanflow

Leading eigenvalue of the linearised Navier-Stokes operator about the mean flow of the shedding Re = 100 cylinder wake.
`../310_stability_direct/direct` solves the same problem about the steady flow. Only the base flow differs.

## Method

- Mode `userParam01 = 3.1` (direct eigensolver, Krylov-Schur), `userParam03 = 2` (two eigenvalues, one complex pair), `userParam07 = 100` (size of the Krylov basis).
  The mesh, `1cyl.usr`, `SIZE` and the settings of the time stepper are those of `../310_stability_direct/direct`. That stage uses `userParam07 = 50` on purpose; the mean flow needs 100.
- `endTime = 1` is the sampling period of the Krylov sequence: one Arnoldi step integrates the linearised equations over one time unit (105 steps of 9.5e-3).
- Base flow: `BF_1cyl0.f00001`, a copy of `../140_baseflow_meanflow/BFA1cyl0.f00001` (`run.local.slurm` copies it). The field files are not tracked.
- No sponge: the lines that set it in `1cyl.usr` are commented out, as in stage 310.

## Result (nd, `ifx`, 8 ranks)

The leading pair has converged with the residual 7.3e-7 after 11 restarts. Wall time 180 s.

| flow | σ | ω | source |
|---|---|---|---|
| steady flow | 0.1248 | 0.7337 | `../310_stability_direct/direct`, this repository |
| mean flow | **0.0022** | **1.0290** | this stage |
| DNS, first harmonic | - | 1.0417 (2π/period of `../140_baseflow_meanflow`) | this repository |
| steady flow | 0.1257 | 0.7389 | Marquet & Carini (2020), Table 2 |
| mean flow | 0.0020 | 1.0322 | Marquet & Carini (2020), Table 2 |
| DNS, first harmonic | - | 1.0446 | Marquet & Carini (2020), Table 1 |

The mean flow is almost marginal: σ is 1.7 % of the growth rate of the steady flow. Its frequency is 1.2 % below the frequency of the DNS, where
the frequency of the steady flow is 30 % below. The ratio of the mean-flow frequency to the DNS frequency is 0.9878 here and 0.9881 in the
published values. The published values come from another mesh and discretisation (47 000 triangles). The mean-flow values differ from them by
0.2e-3 in σ and 0.3 % in ω, the steady-flow values by 0.7 % in σ and in ω.

The leading mode about the mean flow (figure, a) is a wave along the wake. The zero crossings of v on the axis are 2.6 diameters apart downstream of x = 6, a
wavelength of 5.2, and the largest v on the axis is at x = 2.7.

## Reference

`ref/reference.json` checks the leading eigenvalue of `Spectre_NSd.dat` against this run (σ within 5e-4, ω within 2e-3) and against the published
mean-flow values (σ within 2e-3 of 0.0020, ω within 1 % of 1.0322). It checks that the frequency lies within 2 % of the DNS frequency
and that the residual of the leading pair is below 1e-5. A run about the steady flow (σ = 0.1248, ω = 0.7337) fails five of the six checks.
`python3 ../../../scripts/check_against_ref.py .` runs the check.

`plot.py` writes `plot_meanflow_stability.png`: the transverse velocity v of the real part of the leading direct mode, (a) about the mean flow of this stage and
(b) about the steady flow of `../310_stability_direct/direct`. Both panels have the same window, the same colour scale and the same normalisation (maximum of |v| equal to 1, positive sign there). The subtitles give σ and f = ω/2π of the leading pair
(mean flow 0.0022 and 0.1638, steady flow 0.1248 and 0.1168). `ref/meanflow_stability.png` is the figure of the run named in `ref/reference.json`.
`plot_bf.png` (`ref/bf.png`) is the velocity magnitude of the mean flow, drawn the way stage 310 draws its steady base flow.
`plot_spectrum.png` (`ref/spectrum.png`) is the full Ritz spectrum of this stage as sigma against f, with the steady-flow spectrum of stage 310 in gray. Converged values are stars and the legend gives σ and f of the leading one.

## Run

```bash
mks 1cyl
sbatch run.local.slurm      # 8 ranks; run ../140_baseflow_meanflow first
python3 plot.py
python3 ../../../scripts/check_against_ref.py .
```

## Literature

Barkley (2006), Europhys. Lett. 75, 750. Marquet & Carini (2020), arXiv:2008.10002.
