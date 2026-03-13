# Poiseuille RANS — Turbulent channel flow stability

## Physics
2D turbulent channel flow at Re = 100,000 using the k-tau RANS model (m_id = 4). Variable viscosity from the turbulence model, periodic boundary conditions in the streamwise direction, no-slip at top and bottom walls.

Adapted from the Nek5000 NekExamples [RANSChannel](https://github.com/Nek5000/NekExamples/tree/master/RANSChannel) case.

## nekStab Mode
`userParam01 = 3.1` — Direct stability analysis
- `userParam07 = 90` — Krylov subspace dimension
- `variableProperties = yes`, `stressFormulation = yes`
- Scalars: SCALAR01 (tke), SCALAR02 (omega/tau)

### Finite differences for RANS
For RANS stability analysis, finite differences must be used instead of the analytical linearized solver. This is set in `poiseuille_RANS.usr`:

```fortran
iffindiff = .true. ; call bcast(iffindiff,lsize)
```

The linearized RANS equations involve highly nonlinear source terms from the turbulence model; finite differences provide a general-purpose alternative.

## Workflow

### Step 1: Base flow (SFD — Mode 1.1)
Edit `poiseuille_RANS.par`:
```ini
userParam01 = 1.1     # SFD mode
userParam04 = 0.5     # Cutoff frequency
userParam05 = 0.1     # Damping coefficient
```
Output: `BF_poiseuille_RANS0.f00001`

### Step 2: Direct stability (Mode 3.1)
Edit `poiseuille_RANS.par`:
```ini
startFrom = BF_poiseuille_RANS0.f00001
userParam01 = 3.1     # Direct stability analysis
userParam07 = 90      # Krylov subspace dimension
```
Output: eigenmode fields `dRepoiseuille_RANS0.f0000X`, spectrum `Spectre_NSd.dat`

## Prerequisites
- Base flow from Step 1 (or a restart file)

## Run
```bash
mks poiseuille_RANS        # Compile
nekbmpi poiseuille_RANS N  # Run on N MPI ranks
```

## Expected Output
Eigenvalue spectrum of the turbulent channel flow at Re = 100,000.

## Reference
Based on: https://github.com/Nek5000/NekExamples/tree/master/RANSChannel
