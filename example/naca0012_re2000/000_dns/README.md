# naca0012_re2000/000_dns

DNS stage for the NACA 0012 wake at Re=2000.

The companion steady Newton stage gives a stable fixed point, and
`../310_stability_direct/Spectre_NSd_conv.dat` reports the leading pair
`sigma=-0.0665382`, `omega=+/-8.00579`. This stage tests the nonlinear wake
from a finite-amplitude restart and distinguishes ringdown from a sustained
subcritical oscillation.

## Setup

- Mode: DNS (`userParam01 = 0`)
- Mesh/order: copied from the family root (`4368` elements, `lx1=6`,
  `lpmin=8`)
- Reynolds number: `viscosity = -2000.0`
- Sponge: `userParam08=1`, `userParam09=5`, `userParam10=1.7`
- Time step: fixed `dt = 2.0e-4`; `targetCFL = 0.5`
- Output cadence: `writeControl = runTime`, `writeInterval = 10.0`

The initial condition is `resnaca00120.f00010`, copied from the family root.
Its header records `t=9.0`, and it already contains wake structure from the
allowed Re=2000 root DNS series. The run stops at `t=70.0`, giving about 61
convective units after restart, enough to measure ringdown against the stable
linear eigenvalue.

History probe:

1. `(1, 0, 0)` existing wake probe from the family root

## Run

```bash
printf "naca0012\n%s/\n" "$(pwd)" > SESSION.NAME
mks naca0012
sbatch --ntasks=8 --mem-per-cpu=2000 run.local.slurm
```

The Slurm script writes stdout and stderr to `logfile`.

## Latest Result

Jean Zay, 8 ranks, `ifort`, 2026-10-09: wall time 5882 s (98 min). The final checkpoint is `naca00120.f00007` at `t=70.0`.

The finite-amplitude wake rings down monotonically. The late analysis uses the
`v_y(1,0,0)` residual over `30 <= t <= 70`, after the startup transient has passed:

- Hilbert-envelope decay rate: `sigma = -0.063802`
- window-amplitude-ratio decay rate: `sigma = -0.065368`
- linear value of `../310_stability_direct`: `sigma = -0.0665382`
- the DNS decay rates are 4.1 % and 1.8 % less negative than the linear value
- fitted angular frequency: `omega = 8.066107`
- linear value: `omega = 8.00579`
- relative frequency difference: `0.75 %`
- dominant frequency: `f = 1.283761`

The decay rate and the frequency of the nonlinear ringdown agree with the leading eigenvalue of the
Newton base flow.

Run `python3 plot.py` after the Slurm job finishes. It writes:

- `analysis_metrics.dat`
- `signal.png`
- `spectrum.png`
- `snapshot.png`
- `ref/signal.png`
- `ref/spectrum.png`
- `ref/snapshot.png`

`ref/reference.json` holds fixed values of `analysis_metrics.dat` (written by `plot.py` from the probe history)
for the committed run; `plot.py` does not write the reference.

## Verify

```bash
python3 scripts/check_against_ref.py example/naca0012_re2000/000_dns
```
