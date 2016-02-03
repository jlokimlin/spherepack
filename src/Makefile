
include ../make.inc

MODULE_PATH = $(OBJ_DIRECTORY)

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

BIN_PATH = /usr/local/bin

LINKS = \
-lspherepack -L$(EXTERNAL_LIBRARY_PATH) \
-I$(MODULE_PATH)

EXECUTABLE = ./exe

type_threedimensionalvector.o: \
type_threedimensionalvector.f90
	$(FC) $(FFLAGS) -c type_threedimensionalvector.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

type_gaussiangrid.o: \
type_gaussiangrid.f90
	$(FC) $(FFLAGS) -c type_gaussiangrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

type_spherepackworkspace.o: \
type_spherepackworkspace.f90
	$(FC) $(FFLAGS) -c type_spherepackworkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

type_spherepackwrapper.o: \
type_spherepackwrapper.f90
	$(FC) $(FFLAGS) -c type_spherepackwrapper.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

spherepack_wrapper_mod.o: \
type_threedimensionalvector.o \
type_gaussiangrid.o\
type_spherepackworkspace.o\
type_spherepackwrapper.o \
spherepack_wrapper_mod.f90
	$(FC) $(FFLAGS) -c spherepack_wrapper_mod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

main.o: \
spherepack_wrapper_mod.o \
main.f90
	$(FC) $(FFLAGS) -c main.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJECTS_WITH_PREFIX) $(LINKS)

lib: 
	ar rc libspherepack_wrapper.a $(OBJ_DIRECTORY)/*.o
	mv libspherepack_wrapper.a $(LIB_DIRECTORY)

all:
	make clean
	make exec; ./exe
	make lib

clean:
	rm -f $(OBJ_DIRECTORY)/*.o $(OBJ_DIRECTORY)/*.mod $(LIB_DIRECTORY)/*.a
	rm -f $(EXECUTABLE)
