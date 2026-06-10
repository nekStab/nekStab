# cylinder_re100/310_stability_direct/findiff

**Stage**: Stability -- direct, finite-difference Jacobian
**uparam01**: 3.0
**Subroutine**: findiff_stability
**Status**: scaffold only -- no case files yet

**Migration source** (when wiring): (none -- RANS finite-diff path; canonical use is cylinder_re1m)

**Template**: example/_templates/310_stability_direct/findiff/case.par

**Acceptance** (from spine-geometry-list):
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks <case> builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cylinder_re100/310_stability_direct/findiff/ golden artifact captured

