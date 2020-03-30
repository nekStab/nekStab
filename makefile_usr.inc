 $(OBJDIR)/x_krylov.o   :$(CASEDIR)/x_krylov.f;            $(F77) -c $(FL2) -I./ $< -o $@
