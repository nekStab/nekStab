# lid_driven_re3600 / 210_baseflow_newton/fp — Newton steady base flow

**Stage**: 210_baseflow_newton (fp) · **uparam01**: 2 (Newton-GMRES fixed point)
· **Re**: 3600 · **Ranks**: 4. Stage rationale: repo-root `EXAMPLES.md`.

## Seed provenance
`startFrom = BF_cav0.f00001` — a coarse DNS snapshot from the `../../000_dns/`
settle (same `cav` mesh). Re = 3600 is **below** the lid-driven cavity's first
Hopf (~Re 8000), so the flow is steady: plain DNS relaxes toward the fixed
point and Newton-GMRES then converges to it from that native, same-mesh seed
(no foreign field — see the canonical recipe in
`../../../back_fstep_re500/210_baseflow_newton/README.md`).

## Run
```bash
mks cav
sbatch run.local.slurm
```
Sponge disabled (closed/internal geometry). The companion `../upo/` stage is
not applicable here (no periodic orbit at this sub-critical Re).
