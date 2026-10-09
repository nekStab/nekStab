# moving_cylinder_re100/000_dns

DNS stage for the moving-cylinder Re=100 case.

## Setup

- Mode: DNS (`userParam01 = 0`)
- Mesh/order: copied from the family root
- Reynolds number: `viscosity = -100.0`
- Motion/forcing parameters: preserved from the family-root `1cyl.par`
  (`userParam05 = 0.0`, `userParam06 = 0.1643000`)
- Runtime: the case `.usr` sets the final time to `100 / userParam06`
  periods and writes every `10 / userParam06` periods
- History probe: `(5, 0, 0)`, copied from the family-root `1cyl.his`

## Run

```bash
printf "1cyl\n%s/\n" "$(pwd)" > SESSION.NAME
mks 1cyl
sbatch --ntasks=8 --mem-per-cpu=2000 run.local.slurm
```

The Slurm script writes stdout and stderr to `logfile`.

Run `python3 plot.py` after the Slurm job finishes. It writes:

- `signal.png`
- `spectrum.png`
- `snapshot.png`
- `analysis_metrics.dat`
- `ref/signal.png`
- `ref/spectrum.png`
- `ref/snapshot.png`

`ref/reference.json` holds fixed values (the dominant frequency and the late amplitude of the probe signal of the
committed run). `plot.py` reads the probe history and does not write the reference.

## Verify

```bash
python3 scripts/check_against_ref.py example/moving_cylinder_re100/000_dns
```
