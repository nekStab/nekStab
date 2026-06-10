# cubic_cavity/330_transient_growth

**Stage**: transient growth
**uparam01**: 3.3
**Subroutine**: transient growth
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/330_transient_growth/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cubic_cavity/330_transient_growth/ golden artifact captured

