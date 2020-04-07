 $(OBJDIR)/x_eigensolvers.o   :$(NEKSTAB_SOURCE_ROOT)/core/x_eigensolvers.f;      $(F77) -c $(FL2) $< -o $@
 $(OBJDIR)/x_linalg.o         :$(NEKSTAB_SOURCE_ROOT)/core/x_linalg.f;            $(F77) -c $(FL2) $< -o $@
 $(OBJDIR)/x_fixed_point.o    :$(NEKSTAB_SOURCE_ROOT)/core/x_fixed_point.f;       $(F77) -c $(FL2) $< -o $@
 $(OBJDIR)/x_usr_extra.o      :$(NEKSTAB_SOURCE_ROOT)/core/x_usr_extra.f;         $(F77) -c $(FL2) $< -o $@
 $(OBJDIR)/x_utilities.o      :$(NEKSTAB_SOURCE_ROOT)/core/x_utilities.f;         $(F77) -c $(FL2) $< -o $@
 $(OBJDIR)/x_IO.o             :$(NEKSTAB_SOURCE_ROOT)/core/x_IO.f;                $(F77) -c $(FL2) $< -o $@
 $(OBJDIR)/x_postprocessing.o :$(NEKSTAB_SOURCE_ROOT)/core/x_postprocessing.f;    $(F77) -c $(FL2) $< -o $@
