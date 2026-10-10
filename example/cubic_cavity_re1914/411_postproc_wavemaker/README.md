# cubic_cavity_re1914 / 411_postproc_wavemaker

**uparam01**: 4 (post-processing of existing modes: wavemaker 4.1, energy budget 4.2, sensitivity to base-flow changes 4.3)
**Re**: 1914, 16 ranks (`lpmin = 16`). The run takes no time steps and finishes in 3 s.

Structural sensitivity (wavemaker) of the steady cubic cavity flow at Re = 1914, uniform lid. The wavemaker
is the product of the amplitude of the direct mode and the amplitude of the adjoint mode of the leading pair
(sigma = -1.8e-4, omega = 0.587), after the adjoint mode is biorthogonalized to the direct mode (Giannetti
and Luchini, J. Fluid Mech. 581, 2007). It marks where a local change of the flow moves the eigenvalue most.

## Inputs
The run reads these files from its own folder. None of them is tracked: each file is 6 to 12 MB.
- Base flow `BF_cav0.f00001`: copy it from `../210_baseflow_newton/fp/`.
- Direct modes `dRecav0.f0000{1,2}` and `dImcav0.f0000{1,2}`: copy them from `../310_stability_direct/`.
- Adjoint modes `aRecav0.f0000{1,2}` and `aImcav0.f0000{1,2}`: copy them from `../320_stability_adjoint/`.

The code reads the real and the imaginary part of the direct mode from `f00001` and those of the adjoint
mode from `f00002`.

## Run
```bash
mks cav
sbatch run.local.slurm
python3 plot.py
python3 scripts/check_against_ref.py example/cubic_cavity_re1914/411_postproc_wavemaker
```

## Result
Gate run, Jean Zay, 2026-10-10, 16 ranks, with the modes of the gate runs of 310 and 320:
- The wavemaker `wm_cav0.f00001` has its maximum, 45.0, at (x, y, z) = (-0.365, -0.244, 0). A second maximum,
  44.8, lies at (-0.307, 0.333, 0), so the position of the global maximum is not a stable quantity.
- On the plane z = 0 the wavemaker forms a band along the left wall and the bottom wall, with a second
  region under the lid at 0.1 < x < 0.4. The direct mode lies along the left and bottom walls, the adjoint
  mode under the lid and along the left and right walls; the wavemaker is large where the two overlap.
- The run also writes the base-flow sensitivity fields `sr`, `si`, `pr`, `pi`, `tr`, `ti` and the
  energy-budget fields `K01`, `K02` and `PKE_dRecav0.f0000{1,2}`.

## Reference
`ref/reference.json` checks the maximum and the norm of the wavemaker field (1e-3). They are regression
values of this repository; no published value exists for this case. `ref/wavemaker.png` is the figure of the
run named in `ref/reference.json`.

## Figure
`python3 plot.py` writes `plot_wavemaker.png`: (a) the wavemaker, (b) the amplitude of the leading direct
mode and (c) the amplitude of the leading adjoint mode, all on the plane z = 0 and scaled to their maximum.
