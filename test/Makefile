
include ../make.inc

SPHEREPACK_WRAPPER_MODULE_PATH = $(LIB_DIRECTORY)

LINKS = \
-L$(LIB_DIRECTORY) -lspherepack_wrapper \
-I$(SPHEREPACK_WRAPPER_MODULE_PATH) \
-I$(OBJ_DIRECTORY)

EXECUTABLE = ./tests.exe 

#advec.o: \
#advec.f90
#	$(FC) $(FFLAGS) -c advec.f90 $(LINKS)\
#	 -o $(OBJ_DIRECTORY)/$@
#
#helmsph.o: \
#helmsph.f90
#	$(FC) $(FFLAGS) -c helmsph.f90 $(LINKS)\
#	 -o $(OBJ_DIRECTORY)/$@

#shallow.o: \
#shallow.f90
#	$(FC) $(FFLAGS) -c shallow.f90 $(LINKS)\
#	 -o $(OBJ_DIRECTORY)/$@

test.o: \
test.f90
	$(FC) $(FFLAGS) -c test.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@
	
exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

run:
	$(MAKE) clean
	$(MAKE) exec; $(EXECUTABLE)

clean:
	rm -f $(OBJS) $(EXECUTABLE)
