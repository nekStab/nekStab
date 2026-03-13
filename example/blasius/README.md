# Blasius Boundary Layer — Flat plate DNS

## Physics
2D Blasius boundary layer at Re = 740 (based on inlet displacement thickness). DNS of the laminar flat-plate flow.

## nekStab Mode
DNS (legacy .rea format, no nekStab mode parameter).

## Prerequisites
None (self-contained).

## Run
```bash
mks blasius        # Compile
nekbmpi blasius N  # Run on N MPI ranks
```

## Expected Output
Laminar boundary layer profile matching the Blasius similarity solution. Outputs displacement thickness delta* at each time step.

## Reference
Standard boundary layer benchmark (Schlichting, Boundary Layer Theory).
