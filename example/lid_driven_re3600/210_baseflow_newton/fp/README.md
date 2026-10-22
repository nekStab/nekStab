# lid_driven_re3600 / 210_baseflow_newton/fp — Newton steady base flow

**Stage**: 210_baseflow_newton (fp) · **uparam01**: 2 (Newton-GMRES fixed point)
· **Re**: 3600 · **Ranks**: 4. Stage rationale: repo-root `EXAMPLES.md`.

## Seed provenance
`startFrom = BF_cav0.f00001` — the last snapshot (`cav0.f00006`, t = 600) of `../../000_dns/`,
copied to `BF_cav0.f00001` (same `cav` mesh, same order `lx1 = 12`). Re = 3600 is **below** the
lid-driven cavity's first Hopf (~Re 8000), so the flow is steady: plain DNS relaxes toward the
fixed point and Newton-GMRES then converges to it from that native, same-mesh seed (no foreign
field — see the canonical recipe in `../../../back_fstep_re500/210_baseflow_newton/README.md`).
The restart file carries its coordinates, so the cavity height of the seed (`userParam10 = 1.5`)
must be the one of this stage.

The Newton map is the unfiltered time-stepper, while the DNS runs with `filtering = explicit`. At
`lx1 = 12` the seed is a fixed point of the unfiltered map to 4e-9, so Newton needs two
iterations. At lower order the two differ more (table in `../../README.md`).

## Run
```bash
mks cav
sbatch run.local.slurm
```
Sponge disabled (closed/internal geometry). The companion `../upo/` stage is
not applicable here (no periodic orbit at this sub-critical Re).
