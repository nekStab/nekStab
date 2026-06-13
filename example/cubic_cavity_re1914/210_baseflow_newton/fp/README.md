# cubic_cavity_re1914 / 210_baseflow_newton/fp

**Stage**: Newton-Krylov steady baseflow (fixed point)
**uparam01**: 2.0 (steady Newton-GMRES, k_dim = 100)
**Re**: 1914 — just below the first Hopf bifurcation (Re_c ≈ 1916.63;
Loiseau-Robinet-Leriche 2016), so a steady state still exists and is the
physically meaningful base flow here.

**Start from**: `BF_cav0.f00001`, the Re=2500 continuation seed; Newton refines
it down to the Re=1914 fixed point.

**Produces**: the converged `BF_cav0.f00001` steady baseflow — the initial
condition consumed by the downstream stages of this geometry (310 direct,
320 adjoint, 330 transient growth, 411 wavemaker, 600/610/620 modal). Floquet
(311/321) and UPO (210/upo) are not applicable below the Hopf: there is no
periodic orbit to analyse.

**Run**: `sbatch run.local.slurm` (16 ranks). Stage rationale is shared in
`../../../EXAMPLES.md`.

**Acceptance**:
- [x] cav.par populated (uparam01=2.0, viscosity=-1914, k_dim=100)
- [x] cav.usr, cav.re2, cav.ma2, SIZE present (mesh hardlinked, infra copied from root)
- [x] run.local.slurm (16 ranks; writes SESSION.NAME)
- [ ] mks builds clean
- [ ] sbatch run converges to steady baseflow (|F| < 1e-9)
- [ ] refs/ golden artifact captured
