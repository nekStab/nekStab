# Pilot driver — submit_case.py + check_case.py

## What these do

`submit_case.py` runs pre-flight checks on a Nek5000 case directory, stamps
`SESSION.NAME`, submits via `sbatch run.local.slurm`, and records the job in
`validation/campaign_state.json`. `check_case.py` polls Slurm for a submitted
job, inspects new field files on completion, optionally runs plot scripts, and
moves the entry from `active_jobs` to `completed_jobs` or `failed_jobs`.

## Usage

```bash
# Submit (dry run first, then real)
python3 scripts/submit_case.py example/cylinder_re100/210_baseflow_newton/fp --dry-run
python3 scripts/submit_case.py example/cylinder_re100/210_baseflow_newton/fp

# Poll / post-process
python3 scripts/check_case.py  example/cylinder_re100/210_baseflow_newton/fp
python3 scripts/check_case.py  example/cylinder_re100/210_baseflow_newton/fp --plots --strict
```

`--dry-run` skips the actual sbatch call.
`--force` on submit_case.py skips the executable + IC existence checks.
`--strict` on check_case.py exits non-zero on FAILED or AMBIGUOUS outcome.

## State file

`validation/campaign_state.json` is the only mutable state during the pilot.
It holds three lists — `active_jobs`, `completed_jobs`, `failed_jobs`.
Both scripts use `fcntl.flock` for atomic read/write so concurrent submits are safe.

## What is NOT in scope yet

- `validation/catalog.py` status field is not auto-updated; that cross-linking
  is planned for a future phase.
- No automatic re-submission on failure.
- No email/webhook notifications.
