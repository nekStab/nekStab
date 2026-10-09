# Cubic Cavity Re=1950 — unstable periodic limit cycle: Floquet, DNS and modal analysis

Stage rationale: see [EXAMPLES.md](../../EXAMPLES.md).

## Case note

3D lid-driven cubic cavity with a uniform lid, just above the primary Hopf
bifurcation (`Re_c = 1916.6`). The flow settles onto a periodic limit cycle
(period `T = 10.78909`). The cycle is itself unstable at Re = 1950: a pair of
Floquet multipliers lies outside the unit circle. This is the time-periodic
companion to `cubic_cavity_re1914` (the steady-base case that runs the full
SFD / Newton / stability / OTD pipeline). Here the orbit comes directly from
DNS, with no Newton-UPO, and feeds a Floquet stability analysis, a DNS check of
the unstable pair and snapshot-modal decompositions. Thesis source: chapter 4,
`cavsubsec:Non-linearEvolution`, problem definition in
`cav:sec:problem_formulation`.

## Physics
3D lid-driven cubic cavity at Re = 1950 (`viscosity = -1950`), uniform lid
velocity (no regularisation of the edge singularity).

## Workflow (DNS-orbit Floquet, no Newton)
1. **`000_dns/`** — settle onto the limit cycle from a near-periodic start.
2. **`000_dns_period/`** — long DNS + centre probe + FFT for a sharp period
   `T = 10.789`; stamp an on-orbit snapshot as the Floquet base flow.
3. **`311_stability_direct_floquet/`** — direct Floquet (`userParam01 = 3.11`)
   over one period `T`. **`321_stability_adjoint_floquet/`** — the adjoint run.
4. **`350_dns_unstable_orbit/`** — DNS of 2000 time units from the orbit plus
   the Floquet mode, compared with the Floquet pair.
5. **`600_modal_pod/`**, **`610_modal_dmd/`** — POD / DMD of the cycle.

Why DNS-orbit instead of Newton-UPO: the DNS settles directly onto the cycle
(it is unstable only to a slowly growing mode, sigma = 2.9e-3), so Newton
refinement is unnecessary; and a 3D Newton-UPO at a near-marginal Re is stiff
(cf. the `cylinder_re180` 2D-stiffness lesson). The period stamp closes the
loop: see `000_dns_period/README.md`.

## Result — the limit cycle is UNSTABLE to an oscillatory mode
Direct Floquet (`k_dim = 100`, 69 modes converged):
- **Trivial multiplier `mu = 1.0000000`** — the phase mode along the orbit,
  present with residual 2e-16. It confirms `T = 10.78909`.
- **`mu = 0.0121 ± 1.0318i`, `|mu| = 1.0319`** — a pair outside the unit circle:
  `sigma = ln|mu|/T = +2.9e-3`, `omega = 0.1445` (period 43.5 time units, about four
  orbit periods). The adjoint run reproduces it (sigma +2.881e-3, omega 0.14457).
- **The DNS agrees.** Started from the orbit plus 1e-3 of the mode, the
  stroboscopic deviation grows with `sigma = +3.3e-3` and `omega = 0.1423`
  (frequency within 1.6 %, growth rate about 15 % above the Floquet value) and
  saturates near t = 900; the spectrum has the sidebands `n f_1 ± f_F`. Details
  and figure: [`350_dns_unstable_orbit/README.md`](350_dns_unstable_orbit/README.md).

There is a mode to budget, so the energy budget (4.11) is possible on this
cycle. It is exercised on `flip_flop_re62` and `tpjet_re2005`, and is not a
stage here.

POD / DMD: the snapshot stages are described in their own READMEs. The cycle
saturates slowly because Re = 1950 is only just supercritical, and a weak
residual transient can outrank the fundamental in the DMD norm ordering; see
[`610_modal_dmd/README.md`](610_modal_dmd/README.md).

## Run
```bash
mks cav        # compile
sbatch run.local.slurm   # per stage
```
Wall times on 16 ranks: Floquet direct 12.9 h, adjoint 8.4 h, DNS of the unstable
orbit 13.7 h.

## Reference
Standard 3D lid-driven cubic cavity benchmark; first Hopf near Re ≈ 1916.
