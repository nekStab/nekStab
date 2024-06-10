include /home/rfrantz/nekStab/core/makefile_nekStab
$(OBJDIR)/torus_utils.o: $/torus_utils.f90 ;$(FC) -c $(FL2) $< -o $@
