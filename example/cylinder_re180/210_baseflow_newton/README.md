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

- 2D cylinder, Re=180 (`viscosity = -180.0` in `.par`)
- Seed `rstcyl0.f00001`: copy `../000_dns_seed/1cyl0.f00003` (the t = 300
  checkpoint of the Re=180 DNS from rest, on the limit cycle). Field files
  are not tracked; run `000_dns_seed` first (17 min). `startFrom = ... time=0`.
- Period guess `endTime = 5.19022`, from the wake probe of the same DNS
  (`fft_period.py`, 28 cycles).
- `userParam01 = 2.1` — Newton-GMRES UPO shooting
- `userParam07 = 200` — Krylov subspace dimension
- Time integration: auto dt at `targetCFL = 0.5` (`dt = 0`,
  `variableDt = yes`, `timeStepper = bdf3`). The shooting integrates to
  exactly `time = T`; 1190 steps per period.
- `ifdyntol = .false.` in `1cyl.usr`: fixed 1e-9 solver tolerances (see the
  comment there for the run that showed why).

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

## Result

Five Newton iterations from the DNS seed (`residu_newton.dat`):

| Newton iter | residual |
|---|---|
| 1 | 1.2e-2 |
| 2 | 2.8e-3 |
| 3 | 1.5e-5 |
| 4 | 1.2e-9 |
| 5 | 7.8e-10 |

Period **T = 5.189628136**, St = 0.192692034. The DNS period (5.19012 by FFT)
differs by 1e-4 relative. Wall time 54 min on 16 ranks (c8).

`BF_1cyl0.f00001` (written by the run) is the converged orbit. Copy it into
the downstream stages. Its time stamp is the period (nekStab reads the
Floquet period from it). It is the base flow of
`../311_stability_direct_floquet/`, `../321_stability_adjoint_floquet/` and
`../410_postproc_animate_modes/upo/`.

**The downstream stages must use the same `SIZE` as this stage** (polynomial
order 7, `lx1 = 8`, `lxd = 12`). The orbit file holds 8x8 points per
element; at order 5 it is interpolated to a coarser mesh and no longer
closes.

No filter (`filtering = none`): the Floquet stages run without one, so the orbit
comes from the same equations. A comparison with the filter on at order 7 was
not run. (An order-5 comparison moved the multipliers by 1.5 %, but that run
had the mesh mismatch described above.)

Figures: `plot_ic.png` (seed), `plot_residual.png` (GMRES and Newton
residuals), `plot_bf.png` (orbit at t = 0).

## Seeds that failed (history)

- A Re=150 limit cycle with the guess T = 5.5: the period correction
  diverged (T = 538 after one step). The seed was too far from the
  Re=180 orbit.
- Re=170 and Re=175 DNS seeds with their own periods (5.251, 5.221) and
  dynamic solver tolerances: the GMRES solve diverged (3.3e-2 to 9.9) and
  Newton stalled near 1e-2.
- With the Re=180 DNS seed, the Re=180 period and fixed tolerances, Newton
  converges in 5 iterations (table above).

## References

- Barkley, D. & Henderson, R. D. (1996). *Three-dimensional Floquet
  stability analysis of the wake of a circular cylinder.* JFM 322, 215.
- Williamson, C. H. K. (1996). *Vortex dynamics in the cylinder wake.*
  ARFM 28, 477.
- Schuh-Frantz thesis (2022), §appendix:strategies for Newton-UPO
  shooting.
