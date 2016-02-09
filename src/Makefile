
include ../make.inc

MODULE_PATH = $(LIB_DIRECTORY)

LINKS = \
-lspherepack -L$(EXTERNAL_LIBRARY_PATH) \
-I$(MODULE_PATH)

EXECUTABLE = ./exe

type_TrigonometricFunctions.o: \
type_TrigonometricFunctions.f90
	$(FC) $(FFLAGS) -c type_TrigonometricFunctions.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_ThreeDimensionalVector.o: \
type_ThreeDimensionalVector.f90
	$(FC) $(FFLAGS) -c type_ThreeDimensionalVector.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_SphericalUnitVectors.o: \
type_ThreeDimensionalVector.o \
type_SphericalUnitVectors.f90
	$(FC) $(FFLAGS) -c type_SphericalUnitVectors.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_GaussianGrid.o: \
type_GaussianGrid.f90
	$(FC) $(FFLAGS) -c type_GaussianGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_SpherepackWorkspace.o: \
type_SpherepackWorkspace.f90
	$(FC) $(FFLAGS) -c type_SpherepackWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_SpherepackWrapper.o: \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_GaussianGrid.o\
type_SpherepackWorkspace.o\
type_SphericalUnitVectors.o \
type_SpherepackWrapper.f90
	$(FC) $(FFLAGS) -c type_SpherepackWrapper.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

module_spherepack_wrapper.o: \
type_ThreeDimensionalVector.o \
type_GaussianGrid.o\
type_SpherepackWorkspace.o\
type_SpherepackWrapper.o \
module_spherepack_wrapper.f90
	$(FC) $(FFLAGS) -c module_spherepack_wrapper.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
module_spherepack_wrapper.o \
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
