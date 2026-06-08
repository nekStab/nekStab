# Reference results (`ref/`) convention

Every example stage ships a `ref/` folder so a user can **reproduce** the case on
their own machine and **verify** their output against a known-good baseline.

## Layout

```
example/<case>/<stage>/
├─ INPUTS — everything needed to RUN the case
│   *.par  *.usr  SIZE  SESSION.NAME      config
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

| Committed (reproduction ingredients) | Regenerated locally (never committed) |
|---|---|
| config: `*.par *.usr SIZE SESSION.NAME` | `logfile`, `*.log*` |
| mesh: `*.re2 *.ma2` | build artifacts: `nek5000`, `obj/`, `*.o`, `makefile`, `libnek5000.a` |
| IC / seed checkpoints (`startFrom`; `<seed>*.f00001`, `BF_*`, `rst*`) | bulk run output (`<session>0.f0000N`), `residu.dat`, `*.his` |
| `ref/reference.json`, `ref/*.png` | `*.npz` snapshot archives, `*.state` |

Rationale: source lives in `src/` and is compiled per case via `mks`, so build
artifacts are never committed. Field dumps and residual histories are large and
regenerable; their *physics* is captured compactly in `reference.json` + figures.

## `reference.json`

Tolerance-based, not bitwise — spectral-element output is not bit-reproducible
across compiler / MPI / architecture, so verification compares scalars within a
tolerance.

```json
{
  "case": "cylinder_re100/110_baseflow_sfd/akervik",
  "produced": {"nekstab": "v2.0.0-rc2", "nproc": 8, "compiler": "gfortran-13", "date": "2026-06-08"},
  "quantities": {
    "final_residual": {"value": 9.99e-11, "tol_rel": 5e-2, "source": "residu.dat:col2:last"}
  }
}
```

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
