# Cubic Cavity Re=1950 — periodic limit cycle: Floquet + modal analysis

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

3D lid-driven cubic cavity just above the primary Hopf bifurcation
(`Re_c ≈ 1916`), so the flow settles onto a **stable periodic limit cycle**
(period `T ≈ 10.79`). This is the time-periodic companion to
`cubic_cavity_re1914` (the steady-base case that runs the full SFD / Newton /
stability / OTD pipeline). Here the orbit comes directly from DNS — no
Newton-UPO — and feeds a Floquet stability analysis and snapshot-modal
decompositions. Thesis source: chapter 4, `cavsubsec:Non-linearEvolution`,
problem definition in `cav:sec:problem_formulation`.

## Physics
3D lid-driven cubic cavity at Re = 1950 (`viscosity = -1950`). The post-Hopf
state is a stable limit cycle; the DNS converges to it at constant amplitude.

## Workflow (DNS-orbit Floquet — no Newton)
1. **`000_dns/`** — settle onto the limit cycle from a near-periodic start.
2. **`000_dns_period/`** — long DNS + centre probe + FFT for a sharp period
   `T = 10.789`; stamp an on-orbit snapshot as the Floquet base flow.
3. **`311_stability_direct_floquet/`** — direct Floquet (`userParam01 = 3.11`)
   over one period `T`.
4. **`600_modal_pod/`**, **`610_modal_dmd/`** — POD / DMD of the cycle.

Why DNS-orbit instead of Newton-UPO: the DNS settles directly onto the
*stable* cycle, so Newton refinement is unnecessary; and a 3D Newton-UPO at a
near-marginal Re is stiff (cf. the `cylinder_re180` 2D-stiffness lesson). The
period stamp is what closes the loop — see `000_dns_period/README.md`.

## Result — the limit cycle is STABLE
Direct Floquet (k_dim = 192, 170/192 modes converged):
- **Trivial multiplier `μ = 1.000000` recovered to 1e-13** (`σ ≈ 0`) — the
  phase/time-shift mode along the orbit. Its exactness validates `T = 10.789`.
- **Every other `|μ| < 1`** (spectrum decays to `|μ| ≈ 3e-4`): no Floquet mode
  outside the unit circle, consistent with the DNS settling onto the cycle.

There is therefore **no unstable mode to budget**, so this case has no
energy-budget (4.11) stage — that diagnostic is exercised by `flip_flop_re62`
and `tpjet_re2005`, whose cycles are genuinely unstable.

POD: leading real+imag pair holds ≈ 86–94 % of the energy. DMD: eigenvalues lie
near the unit circle and the St = 0.0927 fundamental (matching the Floquet
period) is present — but because Re = 1950 is only just supercritical
(Re_c ≈ 1916) the cycle saturates very slowly, and a weak residual transient
outranks the fundamental in the DMD norm ordering. This is a near-critical-flow
property, not a solver error — see [`610_modal_dmd/README.md`](610_modal_dmd/README.md).

## Run
```bash
mks cav        # compile
sbatch run.local.slurm   # per stage
```

## Reference
Standard 3D lid-driven cubic cavity benchmark; first Hopf near Re ≈ 1916.
