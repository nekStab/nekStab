# cubic_cavity/411_postproc_wavemaker

**Stage**: wavemaker / energy budget
**uparam01**: 4.11
**Subroutine**: stability_energy_budget
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/411_postproc_wavemaker/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cubic_cavity/411_postproc_wavemaker/ golden artifact captured

