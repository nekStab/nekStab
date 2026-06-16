# cubic_cavity/610_modal_dmd

**Stage**: modal DMD
**uparam01**: postproc
**Subroutine**: modal_dmd
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/610_modal_dmd/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cubic_cavity/610_modal_dmd/ golden artifact captured


**PARKED (Re=1914 steady)**: 610_modal_dmd needs a periodic orbit / snapshot series that does not exist below the Hopf. N/A here; applies at re1950 (supercritical).
