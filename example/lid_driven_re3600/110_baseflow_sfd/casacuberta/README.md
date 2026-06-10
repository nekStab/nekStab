# lid_driven/110_baseflow_sfd/casacuberta

**Stage**: baseflow SFD
**uparam01**: 1.1
**Subroutine**: SFD
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/110_baseflow_sfd/casacuberta/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/lid_driven/110_baseflow_sfd/ golden artifact captured

