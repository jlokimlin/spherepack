
include ../make.inc

SPHEREPACK_WRAPPER_MODULE_PATH = $(LIB_DIRECTORY)

LINKS = \
-L$(LIB_DIRECTORY) -lspherepack_wrapper \
-L$(EXTERNAL_LIBRARY_PATH) -lspherepack \
-I$(SPHEREPACK_WRAPPER_MODULE_PATH) \
-I$(OBJ_DIRECTORY)

EXECUTABLE = ./test.exe 

test.o: \
test.f90
	$(FC) $(FFLAGS) -c test.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

run:
	make clean
	make exec; $(EXECUTABLE)

clean:
	rm -f $(OBJS) $(EXECUTABLE)
