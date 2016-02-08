
include ../make.inc

MODULE_PATH = $(LIB_DIRECTORY)

LINKS = \
-lspherepack -L$(EXTERNAL_LIBRARY_PATH) \
-I$(MODULE_PATH)

EXECUTABLE = ./exe

type_trigonometricfunctions.o: \
type_trigonometricfunctions.f90
	$(FC) $(FFLAGS) -c type_trigonometricfunctions.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_threedimensionalvector.o: \
type_threedimensionalvector.f90
	$(FC) $(FFLAGS) -c type_threedimensionalvector.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_sphericalunitvectors.o: \
type_threedimensionalvector.o \
type_sphericalunitvectors.f90
	$(FC) $(FFLAGS) -c type_sphericalunitvectors.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_gaussiangrid.o: \
type_gaussiangrid.f90
	$(FC) $(FFLAGS) -c type_gaussiangrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_spherepackworkspace.o: \
type_spherepackworkspace.f90
	$(FC) $(FFLAGS) -c type_spherepackworkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_spherepackwrapper.o: \
type_trigonometricfunctions.o \
type_threedimensionalvector.o \
type_gaussiangrid.o\
type_spherepackworkspace.o\
type_spherepackwrapper.f90
	$(FC) $(FFLAGS) -c type_spherepackwrapper.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

spherepack_wrapper_mod.o: \
type_threedimensionalvector.o \
type_gaussiangrid.o\
type_spherepackworkspace.o\
type_spherepackwrapper.o \
spherepack_wrapper_mod.f90
	$(FC) $(FFLAGS) -c spherepack_wrapper_mod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

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
	rm -f $(OBJ_DIRECTORY)/*.o $(LIB_DIRECTORY)/*.mod $(LIB_DIRECTORY)/*.a
	rm -f $(EXECUTABLE)
