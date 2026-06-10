# naca0012_re2000/000_dns

DNS stage for the NACA 0012 wake at Re=2000.

The companion steady Newton stage gives a stable fixed point, and
`310_stability_direct/direct/Spectre_NSd.dat` reports the leading pair
`sigma=-0.1127203`, `omega=+/-8.036574`. This stage tests the nonlinear wake
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

Slurm job `184` ran with 8 ranks and completed in `4607 s` of wall time. The
final checkpoint is `naca00120.f00007` at `t=70.0`.

The finite-amplitude wake rings down monotonically. The accepted late analysis
uses the `v_y(1,0,0)` residual over `30 <= t <= 70`, after the startup
transient has passed:

- Hilbert-envelope decay rate: `sigma = -0.063802`
- window-amplitude-ratio decay rate: `sigma = -0.065368`
- linear reference: `sigma = -0.1127203`
- discrepancy: DNS late decay is about `43%` less negative than the Newton
  base-flow eigenvalue
- fitted angular frequency: `omega = 8.066107`
- linear reference: `omega = 8.036574`
- relative frequency error: `0.37%`
- dominant frequency: `f = 1.283761`

The frequency agreement validates the DNS wake period against the linear mode.
The decay-rate discrepancy is not hidden: the DNS steady state differs from the
Newton base flow by O(1e-2) in the wake and O(0.2) in the sponge zone, so the
late DNS ringdown appears to measure the least-damped mode about the DNS steady
state rather than the Newton base-flow operator. This remains an open
investigation.

Run `python3 plot.py` after the Slurm job finishes. It writes:

- `signal.png`
- `spectrum.png`
- `snapshot.png`
- `ref/reference.json`
- `ref/signal.png`
- `ref/spectrum.png`
- `ref/snapshot.png`

The measured outcome and scalar values are recorded in `ref/reference.json`.

## Verify

```bash
python3 scripts/check_against_ref.py example/naca0012_re2000/000_dns
```
