# nekStab Validation Plan

## Overview

Comprehensive validation suite for nekStab stability analysis tools against the AMR paper validation cases. All cases are 2D and test different bifurcation types to demonstrate the solver's capabilities.

**Tolerance**: 5% relative error for quantitative comparisons
**Hardware**: Intel i7-14700 (16 P-cores, pinned via `taskset -c 0-15`)
**Compiler**: Intel ifx 2025.3.2

---

## Validation Cases Summary

| # | Case | Bifurcation Type | Critical Parameter | Status |
|---|------|------------------|-------------------|--------|
| 1 | Cylinder | Hopf | Re_c ≈ 46.6 | ✅ PASSED |
| 2 | Thermosyphon | Pitchfork | Ra_c ≈ 494 | ✅ PASSED |
| 3 | Flip-flop | Neimark-Sacker | Re_c ≈ 61.17 | ❌ INCOMPLETE |
| 4 | Backstep | Transient growth | Re = 500 | ⏳ NOT RUN |
| 5 | Tpjet | Period-doubling | Re_c ≈ 1371 | ⏳ NOT RUN |

---

## Case 1: Cylinder Wake (Hopf Bifurcation)

**Location**: `example/cylinder/stability/direct/`
**Mode**: 3.1 (Direct eigenvalue)
**Test Point**: Re = 50 (above critical Re_c ≈ 46.6)

### Results
```
σ_real = +0.0156    (expected: ~0.016, positive = unstable)
σ_imag = 0.7565     (Strouhal: St = σ_imag/(2π) = 0.1204)
```

### Validation
- [x] Positive growth rate confirms instability above Re_c
- [x] Strouhal number St ≈ 0.12 matches literature
- [x] Complex eigenvalue confirms Hopf bifurcation (oscillatory instability)

**STATUS: ✅ PASSED**

---

## Case 2: Thermosyphon (Pitchfork Bifurcation)

**Location**: `example/thersyphon/stability/direct/`
**Mode**: 3.1 (Direct eigenvalue)
**Test Point**: Ra = 500 (above critical Ra_c ≈ 494)

### Configuration
```ini
[GENERAL]
startFrom = BF_tsyphon0.f00001  # Base flow at Ra=500
userParam01 = 3.1               # Direct eigenvalue
userParam06 = 500.0             # Rayleigh number
userParam07 = 120               # k_dim

[VELOCITY]
viscosity = 5.0                 # Prandtl number Pr=5
```

### Results
```
Leading eigenvalues (Spectre_NSd.dat):
σ₁ = +0.1022 + 0.0i    (UNSTABLE - positive real)
σ₂ = ~0.0 + 0.0i       (marginally stable)
σ₃ = -0.2755 + 0.0i    (stable)
```

### Validation
- [x] Positive real eigenvalue σ = +0.1022 confirms instability
- [x] All eigenvalues are real (σ_imag = 0) → pitchfork bifurcation
- [x] Base flow at Ra=500 correctly shows supercritical instability

**STATUS: ✅ PASSED**

---

## Case 3: Flip-Flop (Neimark-Sacker Bifurcation)

**Location**: `example/flip_flop/baseflow/` (Newton UPO)
**Mode**: 2.1 (Newton-Krylov for UPO)
**Test Point**: Re = 62 (above critical Re_c ≈ 61.17)

### Issue Encountered

Newton-Krylov UPO computation failed with emergency exit:
- Newton iteration 1: Completed (21 GMRES vectors, residual 5.84e-2)
- Newton iteration 2: Started, residual reduced to 2.48e-2, then crashed at step 3258

```
Emergency exit: 3258  time = 8.67765665683918
```

### Analysis

The Newton UPO at Re=62 is challenging because:
1. Starting from Re=60 base flow requires significant orbit correction
2. The period correction was large (~0.87 initially)
3. Numerical instability during Newton iteration 2

### Recommended Solutions

**Option A: Parameter Continuation**
```bash
# Gradual increase: Re=60 → 60.5 → 61 → 61.5 → 62
# Each step uses previous UPO as starting point
```

**Option B: Use Re=60 for Validation (Below Critical)**
```ini
# Run Floquet at Re=60 (below Re_c≈61.17)
# Should show |μ| < 1 (stable multipliers)
userParam01 = 3.11  # Floquet
viscosity = -60.0   # Re=60
startFrom = BF_Re60_2cyl0.f00001
```

**Option C: DNS-based Initial Guess**
```bash
# Run DNS at Re=62, extract approximate UPO from time series
# Use as better initial guess for Newton
```

**STATUS: ❌ INCOMPLETE - Needs parameter continuation or alternative approach**

---

## Case 4: Backward-Facing Step (Transient Growth)

**Location**: `example/back_fstep/transient_growth/`
**Mode**: 3.3 (Transient growth / optimal perturbation)
**Test Point**: Re = 500

### Configuration (to be verified)
```ini
userParam01 = 3.3   # Transient growth mode
userParam07 = 64    # k_dim
```

### Expected Results
- Optimal time τ_opt ≈ 58
- Maximum energy amplification G(τ_opt)

### Files Required
- [x] Base flow: `BF_bfs0.f00001` (exists)
- [ ] Stability config: needs verification

**STATUS: ⏳ NOT YET RUN**

---

## Case 5: Forced Jet (Period-Doubling via Floquet)

**Location**: `example/tpjet/stability/direct_Floquet/`
**Mode**: 3.11 (Direct Floquet)
**Test Point**: Re = 1900 (above critical Re_c ≈ 1371)

### Configuration
```ini
[GENERAL]
startFrom = BF_Re1900_tpjet0.f00001
userParam01 = 3.11    # Direct Floquet
userParam05 = 0.6     # Forcing frequency St=0.6
userParam07 = 128     # k_dim
endTime = 1.666667    # Period = 1/St

[VELOCITY]
viscosity = -1900.0   # Re=1900
```

### Expected Results
- Leading Floquet multiplier μ ≈ -1 (period-doubling)
- |μ| > 1 confirms instability

### Files Required
- [x] Base flow: `BF_Re1900_tpjet0.f00001` (exists)
- [x] Configuration: `tpjet.par` (exists, verified)
- [ ] SIZE optimization for 16 P-cores

**STATUS: ⏳ NOT YET RUN**

---

## Technical Notes

### SIZE File Optimization

For each case, optimize `lpmin` for 16 P-cores:
```fortran
parameter (lpmin=16)              ! 16 P-cores
parameter (lelt=lelg/lpmin + 3)   ! Elements per rank
```

Mesh size constraints:
- Minimum ~20 elements per core for efficiency
- Thermosyphon (lelg=512): use lpmin=8-16
- Flip-flop (lelg=5092): use lpmin=16
- Tpjet: check lelg and adjust

### SESSION.NAME Format
```
casename
/full/path/to/directory/
```

### Running with P-cores Only
```bash
taskset -c 0-15 mpirun -np 16 ./nek5000
```

### Environment Setup
```bash
source ~/.bashrc  # or use 'int' alias
export NEKSTAB_SOURCE_ROOT=/home/rfrantz/nekStab
export NEK_SOURCE_ROOT=/home/rfrantz/nekStab/Nek5000
```

---

## Next Steps

### Priority 1: Complete Remaining Cases
1. [ ] Run backstep transient growth (Re=500)
2. [ ] Run tpjet Floquet (Re=1900)

### Priority 2: Resolve Flip-Flop
3. [ ] Try parameter continuation: Re=60 → 60.5 → 61 → 61.5 → 62
4. [ ] Or validate with Re=60 Floquet (stable case below critical)

### Priority 3: Automation
5. [ ] Create `validate.sh` script for automated testing
6. [ ] Create `validate_compare.py` for result verification
7. [ ] Update `expected_values.json` with actual results

---

## Files Generated

### Thermosyphon Eigenvalues
- `example/thersyphon/stability/direct/Spectre_NSd.dat` - Eigenvalues
- `example/thersyphon/stability/direct/Spectre_NSd_conv.dat` - Convergence

### Flip-Flop (Partial)
- `example/flip_flop/baseflow/residu_newton.dat` - Newton iterations
- `example/flip_flop/baseflow/residu_arnoldi.dat` - GMRES residuals

---

## References

- AMR Paper: Validation cases for nekStab
- Critical values from literature:
  - Cylinder: Re_c ≈ 46.6, St ≈ 0.12
  - Thermosyphon: Ra_c ≈ 494 (pitchfork)
  - Flip-flop: Re_c ≈ 61.17 (Neimark-Sacker)
  - Tpjet: Re_c ≈ 1371 (period-doubling)
  - Backstep: τ_opt ≈ 58 at Re=500
