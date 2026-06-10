# cubic_cavity/620_modal_spod

**Stage**: modal SPOD
**uparam01**: postproc
**Subroutine**: modal_spod
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/620_modal_spod/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cubic_cavity/620_modal_spod/ golden artifact captured

