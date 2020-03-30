 $(OBJDIR)/x_krylov.o   :$(CASEDIR)/x_krylov.f;            $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_stablz.o   :$(CASEDIR)/x_stablz.f;            $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_usrext.o   :$(CASEDIR)/x_usrext.f;            $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_sponge_noise.o   :$(CASEDIR)/x_sponge_noise.f;            $(F77) -c $(FL2) -I./ $< -o $@
 $(OBJDIR)/x_inout.o   :$(CASEDIR)/x_inout.f;            $(F77) -c $(FL2) -I./ $< -o $@
