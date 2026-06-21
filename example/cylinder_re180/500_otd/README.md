# Cylinder OTD — Optimally time-dependent modes

## Physics
2D flow around a circular cylinder at Re = 180. OTD decomposition tracks instantaneous flow instabilities.

## nekStab Mode
`userParam01 = 5` — OTD

## Prerequisites
- Base flow `BF_1cyl0.f00001` (the OTD warm-start, read via `startFrom`). Seed it
  from the Re=180 UPO base flow on the same unified 2128-element mesh:
  `cp ../210_baseflow_newton/BF_1cyl0.f00001 .` — without it the run aborts in
  `mfi_prepare` (missing restart). Field files are gitignored, so this copy is
  required after a fresh checkout.

## Run
```bash
mks 1cyl                     # Compile
cp ../210_baseflow_newton/BF_1cyl0.f00001 .   # seed the OTD base flow
sbatch run.local.slurm       # or: nekbmpi 1cyl N
python3 plot.py              # regenerate the evidence figures
```

## Expected Output
OTD modes and instantaneous Lyapunov exponents (`otd_ftle.dat`,
`otd_growth_rates.dat`, `otd_eigenvalues.dat`) plus mode fields. Krylov
dimension k = 92. Re=180 here is a 2D proxy, so this validates the OTD
machinery, not a physical Mode A/B spectrum.

## Reference
Babaee & Sapsis (2016), Proc. R. Soc. A.
