
include ../compiler_flags.inc
include ./relative_directories.inc

SPHEREPACK_WRAPPER_MODULE_PATH = ../../obj

MODULE_PATH = \
-I$(SPHEREPACK_WRAPPER_MODULE_PATH) \
-I$(OBJ_DIRECTORY)

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

LINKS = \
-L$(LIB_DIRECTORY) -lspherepack_wrapper \
-L$(EXTERNAL_LIBRARY_PATH) -lspherepack \
$(MODULE_PATH)

EXECUTABLE = ./exe

test.o: \
test.f90
	$(FC) $(FFLAGS) -c test.f90 $(LINKS) -o $(OBJ_DIRECTORY)/$@

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJECTS_WITH_PREFIX) $(LINKS)

run:
	make clean
	make exec; ./exe

clean:
	rm -f $(OBJECTS)
	rm -f $(EXECUTABLE)
