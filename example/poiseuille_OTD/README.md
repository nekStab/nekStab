# Poiseuille OTD — Optimally time-dependent modes in channel flow

## Background

This example computes the Optimally Time-Dependent (OTD) decomposition of
2D plane Poiseuille flow at Re = 5000. The OTD method tracks the most
unstable instantaneous directions of the linearized Navier-Stokes operator
by evolving an orthonormal basis that adapts continuously in time.

The nekStab OTD implementation is entirely based on the original code by
**Simon Kern** ([OTD4NekStab](https://github.com/Simkern/OTD4NekStab/tree/main/src_base)).
All credits for the OTD method implementation, algorithm design, and
mathematical formulation go to Simon Kern. The adaptation into nekStab
integrates his work into the perturbation infrastructure, using Nek5000's
native `incompLinNS` equation type and mass-weighted inner products.

### Interpreting the Lyapunov exponents

The OTD decomposition produces **finite-time Lyapunov exponents (FTLEs)**,
which measure the average exponential growth or decay rate of perturbations
along each OTD direction over a time window `[0, T]`:

```
lambda_i(T) = (1/T) integral_0^T  L_r(i,i)(t) dt
```

where `L_r = U^T L_{NS} U` is the reduced operator projected onto the OTD
basis.

- **lambda > 0**: perturbations grow exponentially on average along this
  direction — the flow is instantaneously unstable.
- **lambda < 0**: perturbations decay — the flow is stable along this
  direction.
- **lambda ~ 0**: marginal/neutral — perturbations neither grow nor decay.

For this Poiseuille flow at Re = 5000 (below the critical Re ~ 5772 for
Tollmien-Schlichting instability), both exponents converge to negative
values (lambda_1 ~ -0.012, lambda_2 ~ -0.025), confirming that the flow
is linearly stable: all perturbations decay.

The three output files report different spectral quantities of the reduced
operator `L_r`:

| File | Content |
|------|---------|
| `otd_ftle.dat` | FTLEs: time-averaged diagonal of `L_r` (Lyapunov exponents) |
| `otd_eigenvalues.dat` | Eigenvalues of the full (non-symmetric) `L_r` — instantaneous growth rates |
| `otd_growth_rates.dat` | Eigenvalues of `(L_r + L_r^T)/2` — symmetric part of `L_r` |

The eigenvalues of `L_r` give **instantaneous** growth rates (which fluctuate
in time), while the FTLEs give the **time-averaged** rates (which converge
as T grows). The symmetric part `L_s` isolates the contribution from
strain-like mechanisms (excluding rotation).

## nekStab Mode

OTD via Nek5000's `incompLinNS` equation type (not nekStab modes).
- `numberOfPerturbations = 2` — two OTD basis vectors
- `solveBaseflow = no` — analytical parabolic profile

Key OTD parameters (set in `.usr`):
- `otd_iostep = 50` — output frequency
- `otd_gsstep = 2` — Gram-Schmidt reorthonormalization frequency
- `otd_prntsp = 20` — print frequency

## Prerequisites

None (base flow is the analytical Poiseuille profile `U = 1 - y^2`).

## Run

```bash
mks poiseuille_OTD        # Compile
nekbmpi poiseuille_OTD N  # Run on N MPI ranks
```

## Expected Output

OTD modes and Lyapunov exponents for Poiseuille flow. Both FTLEs converge
to negative values, consistent with the flow being linearly stable at
Re = 5000 (below Re_c ~ 5772).

## Credits

The OTD module in nekStab was developed by **Simon Kern**. The original
implementation is available at https://github.com/Simkern/OTD4NekStab.

## Reference

- Kern, OTD4NekStab: https://github.com/Simkern/OTD4NekStab
- Babaee & Sapsis (2016), *A minimization principle for the description of
  modes associated with finite-time instabilities*, Proc. R. Soc. A 472.
- Babaee, Farazmand, Haller & Sapsis (2017), *A reduced-order description
  of the dynamics in and of FTLEs*, Chaos 27.
- Blanchard & Sapsis (2019), *Analytical description of optimally
  time-dependent modes of decay in the Stokes limit*, Phys. Fluids 31.
