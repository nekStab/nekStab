# cubic_cavity/320_stability_adjoint

**Stage**: stability adjoint
**uparam01**: 3.2
**Subroutine**: Krylov-Schur (adjoint)
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/320_stability_adjoint/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/cubic_cavity/320_stability_adjoint/ golden artifact captured

