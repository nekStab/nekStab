# Slot Jet with Free-Stream Turbulence — DNS

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
