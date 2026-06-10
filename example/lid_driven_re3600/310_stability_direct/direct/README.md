# lid_driven/310_stability_direct/direct

**Stage**: stability direct
**uparam01**: 3.1
**Subroutine**: Krylov-Schur
**Status**: scaffold only — no case files yet

**Migration source**: (none — to build from template)

**Template**: example/_templates/310_stability_direct/direct/case.par

**Acceptance**:
- [ ] case.par populated from template
- [ ] *.usr, *.re2, *.ma2, SIZE present (symlink or copy from migration source)
- [ ] SESSION.NAME set
- [ ] run.local.slurm with 8 ranks
- [ ] mks builds clean
- [ ] sbatch run completes, produces expected artifacts
- [ ] refs/lid_driven/310_stability_direct/ golden artifact captured

