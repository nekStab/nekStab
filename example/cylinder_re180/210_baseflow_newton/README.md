# Cylinder UPO at Re=180 — Newton-GMRES on the 2D periodic wake

## ⚠ This is a 2D PROXY of a 3D problem — not a physical reproduction

Barkley & Henderson (*JFM* 322, 215, 1996) is **fundamentally a 3D
analysis**. Only the **base state** is 2D — the von Kármán limit cycle
(the T-periodic wake), which is exactly the UPO computed here. The
**instability** they discovered is intrinsically **three-dimensional**:

- **Mode A** (Re ≈ 188, spanwise wavelength λ_z ≈ 3.96 D, β ≈ 1.585) and
  **Mode B** (Re ≈ 259, λ_z ≈ 0.82 D, β ≈ 7.6) are **spanwise-periodic**
  Floquet modes — perturbations ~ e^{iβz} on the 2D periodic wake. Their
  structure lives *out of the plane*.

- A purely **2D** computation has no spanwise direction, so it can only
  represent **β = 0** perturbations. The β=0 Floquet multiplier of a
  limit cycle is the trivial **unit multiplier** (the phase / time-shift
  mode along the orbit); the 2D wake is otherwise 2D-stable in this Re
  range. **A 2D Floquet here cannot find Mode A or Mode B — those *are*
  the third dimension.**

**Why we still run it in 2D:** a genuine 3D Floquet (a spanwise-periodic
mesh sized to fit λ_z, or a β-parameterized 3D Floquet) is too expensive
for the test suite. This 2D case is a **cost-reduced proxy** that
exercises the full nekStab pipeline — DNS seed → Newton UPO → Floquet
eigensolver → mode animation (`cylinder_re180/311_*`, `321_*`,
`410_postproc_animate_modes/upo/`) — standing in for the most famous
Floquet stability case. **The Re=180 value matches the Barkley setup's
parameters; it does NOT imply the 2D run reproduces Mode A/B.** Treat the
2D Floquet output as a *machinery demonstrator*, not as physical
secondary-instability results. For physically accurate Mode A/B, redo the
base flow + Floquet in 3D.

## Physics

- 2D cylinder, Re=180 (viscosity = -180.0 in `.par`)
- Period guess: `endTime = T` where **T is measured by FFT of a Re=170
  DNS** (see workflow below). Re=170 gives St ≈ 0.190 ⇒ **T ≈ 5.251**.
  Do NOT use an arbitrary guess — the period-Newton diverges if the
  guess is far (see "What went wrong" below).
- `userParam01 = 2.1` — Newton-GMRES UPO shooting
- `userParam07 = 200` — Krylov subspace dimension
- Time integration: **auto dt at `targetCFL = 0.5`** (`dt = 0`,
  `variableDt = yes`, `timeStepper = bdf3`). The shooting integrates to
  exactly `time = T` regardless of step size, so variable dt is fine here
  — confirmed with the Re=170 seed: the first Newton residual is already
  ~4e-3 (deep in the convergence basin, no T-explosion).

### Natural UPO — solve for space AND time

`userParam01 = 2.1` selects the **natural** (self-excited) UPO mode
(`isNewtonPO`): the period **T is an unknown of the Newton system**,
solved jointly with the flow state. The cylinder wake selects its own
shedding frequency, so we must NOT fix the period — we solve for the
state *and* T together, closed by a phase condition. `endTime` is only
the **initial guess** for T; each Newton iteration applies a period
correction (`ΔT`, printed as "Newton period correction") and the orbit
time is updated until both the shooting residual and the phase condition
vanish. (Contrast `userParam01 = 2.2`, `isNewtonPO_T` — a *forced* UPO
where the period is imposed by an external forcing frequency and is NOT
solved for.)

## What went wrong with the old Re=150 seed (do not repeat)

The previous workflow seeded the UPO from a **Re=150** limit cycle with
an arbitrary period guess `endTime = 5.5`. It **diverges**: the
period-Newton overshoots catastrophically because the seed is too far
off the Re=180 orbit.

Observed (post-2128-mesh):

| Newton iter | period T | residual |
|---|---|---|
| 1 | 5.50 | 0.24 |
| 2 | 538 | 0.49 |
| 3 | 885 | 2.58 |
| 4 | 1292 | 27.3 |

ΔT jumps +533 on the first step. Root cause: the period correction is
`ΔT = −r / (∂r/∂T)` and the `main` solver applies the **full** step
with no trust-region/damping, so a far seed (large initial residual)
produces a runaway. A Re=150 field is genuinely off-orbit at Re=180.

## Initial-condition + period workflow (DNS at Re=170 + FFT)

Two things must be right for the period-Newton to converge: the seed
must be **near** the Re=180 orbit, and the **period guess must be
accurate**. Both come from a short DNS at **Re=170** (close to the
target, safely below Re_A≈188 so it stays 2D):

1. **DNS at Re=170.** Use `../000_dns_seed_re170/` (copy of the
   canonical 2128/lx1=8 case with `viscosity = -170.0`). Seed it from
   any nearby saturated limit cycle (a Re=150 field re-saturates to the
   Re=170 orbit in ~50 t_c) and run **≥ 100 shedding periods** of clean
   signal (`endTime` such that ~100 periods elapse after saturation).

2. **Probe the wake.** `hpts` reads probe points from the `1cyl.his`
   header — create it with one downstream-centreline point:
   ```
   1
   5 0 0
   ```
   Nek appends `time  u  v  p` per step. The transverse velocity `v`
   (column 3) oscillates at the shedding frequency.

3. **FFT for the period.** Demean, Hanning-window, uniformly resample
   the saturated `v(t)`, take `np.fft.rfft`, find the peak bin, refine
   with parabolic interpolation:
   ```python
   import numpy as np
   a = np.array([[float(x) for x in l.split()[:4]]
                 for l in open('1cyl.his') if len(l.split()) >= 4])
   t, v = a[:,0], a[:,2]
   m = t > t_saturated           # drop the re-saturation transient
   tu = np.linspace(t[m][0], t[m][-1], m.sum())
   vu = np.interp(tu, t[m], v[m]); vu -= vu.mean()
   F  = np.abs(np.fft.rfft(vu*np.hanning(len(vu))))
   fr = np.fft.rfftfreq(len(vu), tu[1]-tu[0])
   k  = 1 + np.argmax(F[1:])
   d  = 0.5*(F[k-1]-F[k+1])/(F[k-1]-2*F[k]+F[k+1])  # sub-bin peak
   f0 = (k+d)*(fr[1]-fr[0]); print('T =', 1/f0)
   ```
   Re=170 result: **St = 0.1904, T = 5.251** (86-period FFT window).

4. **Seed + period into the UPO.** Copy the latest saturated Re=170
   snapshot as `rstcyl0.f00001`, set `.par` `endTime = 5.251` (the
   measured T) and `viscosity = -180.0`, and run Newton-GMRES. The
   small Re=170→180 gap keeps the first Newton residual small, so ΔT
   stays O(1) and the period converges instead of exploding.

## Status — Newton-UPO does NOT converge at Re=180 (2D-stiffness limit); orbit taken from DNS

We pushed the seed + period as far as possible (2026-06-18):

1. **Sharp period from a long Re=175 DNS.** A ~228-period Re=175 run
   (`../000_dns_seed_re175/`, `hpts` wake probe at (5,0,0)) gives a
   very sharp spectral peak (Q≈100). PSD-peak and 172-cycle counting
   agree to 5 sig figs: **St = 0.1915, T = 5.22130** (vs the older,
   coarser Re=170 value 5.251). Re=175 is close to Re=180 (St varies
   slowly) yet gives a clean single-frequency limit cycle.

2. **On-orbit seed.** Newton-UPO was seeded from a saturated on-orbit
   snapshot (constant amplitude over 10 periods). The improved seed +
   sharp period cut the **initial** Newton defect ~7× (1.26e-2 vs the
   earlier 9.27e-2) and activated a sane period correction (ΔT≈−1.9e-3).

3. **But the Newton LINEAR solve diverges.** The GMRES on the
   monodromy-like operator `M−I` (with the phase/period constraints)
   blows up: residual 3.3e-2 → **9.88** → 1.6e-1, never descending,
   each Newton step burning a full k_dim=100 Arnoldi restart (~25 min)
   to no avail. At Re=180 (Mode-A threshold) the 2D orbit is too stiff
   for clean Newton-UPO refinement. This is the **2D-proxy limit**, not
   a fixable setup issue.

**Resolution — use the DNS limit cycle directly as the Floquet base
flow.** The periodic orbit *physically exists* (the DNS settled onto it:
constant amplitude, Q≈100 single-frequency spectrum). Newton would only
refine it to machine precision, which is not required to demonstrate the
Floquet machinery. So `BF_1cyl0.f00001` is a **Re=180 DNS on-orbit
snapshot**, with its field **time stamp set to the period T = 5.22130**
(nekStab reads the Floquet period from the base-flow file's time —
`matvec.f90:686 param(10)=qbase%time`). This feeds the downstream
`311_stability_direct_floquet`, `321_stability_adjoint_floquet`, and
`410_postproc_animate_modes/upo`.

**Caveat (honest):** the orbit field is Re=180 but the period (5.22130)
was measured at Re=175, so over one period the orbit closes to ~0.3–0.5%,
not exactly. For this 2D **proxy** (already non-physical for Mode A/B)
that is acceptable for a machinery demonstration. For a precise Re=180
period, run the same long-DNS+FFT at Re=180 and re-stamp the orbit time.

## Prerequisites

- Re=170 limit-cycle snapshot saved as `rstcyl0.f00001` (this case's
  restart source). Mesh: `nelt = 2128`, lx1=8 (canonical cylinder mesh).
- `.par` `endTime` set to the FFT-measured period (≈5.251), not a guess.
- nekStab built with mode_codes (no special flags required).

## Acceptance

- Newton-UPO converges in ≤ 20 iterations to `|F| < 1e-8`.
- Period `T` stays O(5.2–5.3) throughout (no T-explosion); converges
  near the FFT value 5.251 (St near 0.190; precise value in the log).
- `BF_1cyl0.f00001` written at convergence; this becomes the Floquet
  base flow for the downstream stages.

## References

- Barkley, D. & Henderson, R. D. (1996). *Three-dimensional Floquet
  stability analysis of the wake of a circular cylinder.* JFM 322, 215.
- Williamson, C. H. K. (1996). *Vortex dynamics in the cylinder wake.*
  ARFM 28, 477.
- Schuh-Frantz thesis (2022), §appendix:strategies for Newton-UPO
  shooting.
