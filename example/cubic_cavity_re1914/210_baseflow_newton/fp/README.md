# cubic_cavity/210_baseflow_newton/fp

**Stage**: baseflow Newton
**uparam01**: 2.0/2.1
**Subroutine**: newton_krylov
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/210_baseflow_newton/fp/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cubic_cavity/210_baseflow_newton/ golden artifact captured

