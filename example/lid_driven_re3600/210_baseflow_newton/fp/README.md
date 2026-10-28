# lid_driven_re3600 / 210_baseflow_newton/fp — Newton steady base flow

**Stage**: 210_baseflow_newton (fp) · **uparam01**: 2 (Newton-GMRES fixed point)
· **Re**: 3600 · **Ranks**: 4. Stage rationale: repo-root `EXAMPLES.md`.

## Seed provenance
`startFrom = BF_cav0.f00001` — the first snapshot (`cav0.f00001`, t = 100) of `../../000_dns/`,
copied to `BF_cav0.f00001` (same `cav` mesh, same order `lx1 = 12`). Re = 3600 is **below** the
lid-driven cavity's first Hopf (~Re 8000), so the flow is steady: plain DNS relaxes toward the
fixed point and Newton-GMRES then converges to it from that native, same-mesh seed (no foreign
field — see the canonical recipe in `../../../back_fstep_re500/210_baseflow_newton/README.md`).
The restart file carries its coordinates, so the cavity height of the seed (`userParam10 = 1.5`)
must be the one of this stage.

The t = 100 state is a fixed point to 4.3e-7, so Newton has work to do: 10 iterations and 117
Krylov vectors to the tolerance 1e-12 (2026-10-28, Jean Zay, 4 ranks, 136 s). The t = 600 state
is a fixed point to 4.4e-9 and needs 6 iterations to the same tolerance. Both seeds end on the
same field to 8e-5: the 3e-3 gap of the earlier runs came from the tolerance, not from a
second steady state.

The tolerance (`residualTol`) applies to the squared L2 norm of the residual over the unit Newton
time. The slowest mode decays at rate 0.0168, so the field error is about sqrt(tol)/0.0167. At
the tolerance 1e-9, Newton from t = 100 and from t = 200 stopped 3e-3 from the steady state;
the tolerance 1e-12 brings that below 1e-4. `writeDoublePrecision = yes` keeps the base flow
at that accuracy for the stability stages.

The Newton map is the unfiltered time-stepper, while the DNS runs with `filtering = explicit`.
At lower order the two differ more (table in `../../README.md`).

## Run
```bash
mks cav
sbatch run.local.slurm
```
Sponge disabled (closed/internal geometry). The companion `../upo/` stage is
not applicable here (no periodic orbit at this sub-critical Re).
