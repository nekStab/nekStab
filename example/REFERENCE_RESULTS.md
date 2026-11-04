# Reference results (`ref/`) convention

Every validated example stage ships a `ref/` folder so a user can **reproduce** the case on
their own machine and **verify** their output against a known-good baseline. The experimental
`cylinder_re1m` and the smoke test `cylinder_re100/000_dns/ci_test` have no `ref/` folder.

## Layout

```
example/<case>/<stage>/
├─ INPUTS — everything needed to RUN the case
│   *.par  *.usr  SIZE                    config
│   *.re2  *.ma2                          mesh
│   <seed>0.f00001                        initial condition / startFrom checkpoint
│   plot.py  README.md  run.local.slurm
│
└─ ref/  — EXPECTED OUTPUTS, compact and version-controlled
    ├─ reference.json   key scalars + tolerances + provenance
    └─ *.png            reference figures (also read by the validation gallery)
```

The case directory holds **inputs**; `ref/` holds **expected outputs**. The `ref/`
figures are the single source of truth for both user verification and the gallery.

## What is committed vs regenerated

| Committed (reproduction ingredients) | Regenerated locally (not committed) |
|---|---|
| config: `*.par *.usr SIZE` | `SESSION.NAME` (gitignored; `nekbmpi` and `run.local.slurm` write it), `logfile`, `*.log*` |
| mesh: `*.re2 *.ma2` | build artifacts: `nek5000`, `obj/`, `*.o`, `makefile`, `libnek5000.a` |
| IC / seed checkpoints (`startFrom`; `<seed>*.f00001`, `BF_*`, `rst*`) | bulk run output (`<session>0.f0000N`), except the result fields of "What a stage ships" |
| `ref/reference.json`, `ref/*.png` | `*.npz` snapshot archives, `*.state` |
| `.his` and `.dat` files (`residu.dat`) up to 1 MB, as listed in `shipped_results.tsv` | other `.his` and `.dat` files |

Rationale: source lives in `src/` and is compiled per case via `mks`, so build
artifacts are never committed. Field dumps are large and
regenerable; their *physics* is captured compactly in `reference.json` + figures.

## `reference.json`

Tolerance-based, not bitwise — spectral-element output is not bit-reproducible
across compiler / MPI / architecture, so verification compares scalars within a
tolerance.

```json
{
  "case": "cylinder_re100/110_baseflow_sfd/akervik",
  "summary": "Steady base flow of the cylinder wake at Re = 100 by selective frequency damping.",
  "description": "What the quantities check, and where the values come from.",
  "produced": {"nekstab": "v2.0.0-rc2", "nproc": 8, "compiler": "gfortran-13", "date": "2026-06-08"},
  "quantities": {
    "final_residual": {"value": 9.99e-11, "tol_rel": 5e-2, "source": "residu.dat:col2:last"}
  }
}
```

`summary` is the one line on the card of the stage in the validation gallery: what the stage
computes. `description` says what the quantities check. Both are required, and so is at least one
figure in `ref/` and a `README.md` in the stage folder: `validation/tests/test_gallery_complete.py`
fails for a stage folder without them. Run `python validation/build_gallery.py` after a change.

Each quantity names its `source` (how to extract it from a fresh run), an expected
`value`, and a tolerance (`tol_rel` relative or `tol_abs` absolute). Typical
quantities: `final_residual` (baseflow), `leading_eig_real`/`leading_eig_imag`
(stability), `period_T` (UPO), `strouhal`/`cd_mean` (DNS).

## Verifying a run

```
scripts/check_against_ref.py example/<case>/<stage>
```

Re-extracts each `reference.json` quantity from the freshly produced output,
prints `PASS`/`FAIL` per quantity with relative error vs tolerance, and exits
non-zero if any quantity is out of tolerance — usable both as a user smoke test
and a CI gate.

## What a stage ships

These are the rules for the files next to `ref/`, so that the repository stays small and each stage still shows its result.
`example/shipped_results.tsv` lists each file, its md5 and where it came from.

| Stage | Result field |
|:---|:---|
| 000, 010 (DNS) | last snapshot |
| 110 to 140, 210 (base flow) | `BF_*` (140: the mean `BFA*`) |
| 310, 311, 315, 320, 321 (eigenmodes) | real and imaginary part of `f00001`, the leading mode |
| 330 (transient growth) | first optimal: `pRe*`, `pIm*`, `ore*` |
| 410 (animation) | last frame |
| 411, 412 (sensitivity) | `wm_*` or `fsr*`; the budget stages: `F01*` and the `PKE_*` tables |
| 500 (OTD) | `BF_*` and the last `r01*` |
| 600, 610, 620 (POD, DMD, SPOD) | mean and leading mode: first `pod*`, shedding `dm1*` and `dm2*`, saved `sRe*` and `sIm*` of the largest eigenvalue |
| 630 (DMT bank) | last snapshot and filter 1 |

- **Start fields.** A stage also ships the fields it opens at start. The cubic cavity cases do not: their READMEs name the stage to run first and the copy. Git stores identical files once, so a copy in a later stage adds nothing to a clone.
- **Cubic cavity.** The 3D fields are large: text, figures and one start field per case.
- **Text and figures.** All figures. `.his` and `.dat` files up to 1 MB. `residu_dmt.dat` stays out, because a run appends to an existing file.
- **Never:** `logfile`, `.nek5000`, `SPG*`, `KRY*`, snapshot series. `.gitignore` blocks most shipped names, so add a shipped file with `git add -f`.

**File format.** Results and frozen base flows are single precision without pressure. Modes and other results have no coordinates; base flows and other start fields keep them, because `plot.py` and field viewers read them from the file.
The start field of a stage that marches in time (000, 010, 110 to 140, 210, 600 to 630) stays as the run wrote it. A restart without pressure moved `snapshot_vy_max` of `cylinder_re100/000_dns` by 1.1e-6, above its tolerance of 1e-6. With the original field the change was 4e-10. A frozen base flow in this format moved the leading eigenvalues by 1e-7 or less.

**After a run.** A run overwrites these files, and `git status` lists what changed. Fields and figures are binary, so compare the values with `scripts/check_against_ref.py`. It compares only the quantities in `ref/reference.json`, not the other fields or the figures. It refuses a file that still has the shipped bytes, because no run wrote it.

**`plot.py` without a run.** A stage ships only part of what its `plot.py` reads. On the shipped files it redraws what it can and stops at the first missing file, whose name it prints; the cubic cavity stages need a run. A figure redrawn this way can differ from the shipped one, because a mode has an arbitrary sign and the shipped figure can come from an earlier file of the series. Coordinates of a field that has none come from a field of the same folder.
