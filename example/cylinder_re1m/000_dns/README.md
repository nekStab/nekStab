# cylinder_re1m/000_dns

**Stage**: DNS
**uparam01**: 0
**Subroutine**: no-op
**Status**: scaffold only — no case files yet

**Migration source**: example/cylinder/RANS

**Template**: example/_templates/000_dns/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cylinder_re1m/000_dns/ golden artifact captured

