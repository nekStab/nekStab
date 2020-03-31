 $(OBJDIR)/x_eigensolvers.o   :$(CASEDIR)/x_eigensolvers.f;      $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_linalg.o         :$(CASEDIR)/x_linalg.f;            $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_fixed_point.o    :$(CASEDIR)/x_fixed_point.f;       $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_usr_extra.o      :$(CASEDIR)/x_usr_extra.f;         $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_utilities.o      :$(CASEDIR)/x_utilities.f;         $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_IO.o             :$(CASEDIR)/x_IO.f;                $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_postprocessing.o :$(CASEDIR)/x_postprocessing.f;    $(F77) -c $(FL2) -I./ $< -o $@
