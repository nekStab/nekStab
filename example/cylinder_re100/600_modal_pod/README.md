# Cylinder Modal Analysis — POD/DMD/SPOD from DNS snapshots

## Physics
2D flow around a circular cylinder at Re = 100. Modal decomposition of DNS snapshots.

## nekStab Mode
`userParam01 = 6` — Modal analysis (runs POD, DMD, and SPOD together)

Use `6.1`, `6.2`, or `6.3` to run POD, DMD, or SPOD individually.

## Two-Phase Workflow
`run.local.slurm` runs both phases in one job and leaves `1cyl.par` as the modal-analysis file:
1. **Phase 1, snapshots (`dns.par`, `userParam01 = 0`):** DNS from the restart file `rst_1cyl0.f00001` (limit-cycle state at t = 150) to t = 200. It writes 100 field snapshots, one every 0.5 time units (about 100 s on 8 ranks).
2. **Phase 2, modal analysis (`1cyl.par`, `userParam01 = 6`):** reads the snapshots and computes POD, DMD and SPOD. The `.usr` file selects the three methods only when `userParam01 >= 6`.

## Prerequisites
- Restart file: `rst_1cyl0.f00001` (tracked; t = 150).

## Run
```bash
mks 1cyl                  # compile
sbatch run.local.slurm    # both phases, 8 ranks
```

## Result
Jean Zay, 8 and 16 ranks, 2026-10-09 (identical to six digits at both rank counts):
- POD: the leading pair holds 49.87 % and 47.72 % of the fluctuation energy
  (97.59 % cumulative); `pod_spectrum.dat`.
- DMD (rank 20): the shedding pair has St = 0.16566, its harmonics 0.33158 and
  0.66301; `dmd_spectrum.dat`. The 2D cylinder at Re = 100 sheds at St about 0.165.
- Streaming SPOD (nfft = 64, noverlap = 48): the leading eigenvalue peaks in the
  bin St = 0.15625 (15.77) and the next bin (0.1875) holds 8.41; the shedding
  frequency lies between the two bins. In every bin the three eigenvalues sum to the
  total power of the POD-FFT analysis (`pod_fft_spectrum.dat`) to four digits.

`ref/reference.json` holds these values. Check with
`uv run --with numpy --with pymech python scripts/check_against_ref.py example/cylinder_re100/600_modal_pod`.

## Reference
Barkley & Henderson (1996), JFM 322.
