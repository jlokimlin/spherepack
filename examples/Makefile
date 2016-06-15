
include ../make.inc

LINKS = \
-L$(LIB_DIRECTORY) -l$(LIB_NAME) \
-I$(MOD_DIRECTORY) \
-I$(OBJ_DIRECTORY)

advec:
	$(FC) $(FCFLAGS) -o advec.exe type_AdvectionSolver.f90 advec.f90 $(LINKS)
	./advec.exe	

helmsph:
	$(FC) $(FCFLAGS) -o helmsph.exe helmsph.f90 $(LINKS)
	./helmsph.exe	

shallow:
	$(FC) $(FCFLAGS) -o shallow.exe type_ShallowWaterSolver.f90 shallow.f90 $(LINKS)
	./shallow.exe	

run: 
	$(MAKE) clean
	$(MAKE) examples

examples: advec helmsph shallow

clean:
	rm -f $(OBJS) *.exe *.mod
