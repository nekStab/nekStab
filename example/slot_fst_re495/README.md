# Slot Jet with Free-Stream Turbulence — DNS

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

This case is a DNS/inflow demonstration for free-stream-turbulence synthesis in
a slot-jet geometry. The named thesis chapters discuss jet-transition problems,
but they do not contain a dedicated slot-jet FST section. Source for the
methodological framing: thesis chapter 6, `jet_sec:intro` and
`jet_sec:nonlinear`; FST provenance is documented in the project
acknowledgments.

## Physics
2D slot jet at Re = 495 with free-stream turbulence. DNS with OIFS time stepping (CFL = 3) and high-pass filter (hpfrt).

## nekStab Mode
`userParam01 = 0` — DNS
- `userParam06 = 0.35` — jet amplitude R
- `userParam07 = 0.3603` — delta*/D

## Prerequisites
- Restart file: `34_slot0.f00001`

## Run
```bash
mks slot        # Compile
nekbmpi slot N  # Run on N MPI ranks
```

## Expected Output
DNS of the slot jet mixing layer. Vortex roll-up and pairing dynamics.

## Reference
Internal example.
