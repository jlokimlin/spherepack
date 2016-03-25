
include ../make.inc

MODULE_PATH = $(LIB_DIRECTORY)

LINKS = \
-lspherepack -L$(EXTERNAL_LIBRARY_PATH) \
-I$(MODULE_PATH)

EXECUTABLE = ./build_lib.exe

type_SphericalGrid.o: \
type_SphericalGrid.f90
	$(FC) $(FFLAGS) -c type_SphericalGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)
	 
type_RegularGrid.o: \
type_SphericalGrid.o\
type_RegularGrid.f90
	$(FC) $(FFLAGS) -c type_RegularGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)
	 
type_GaussianGrid.o: \
type_SphericalGrid.o\
type_GaussianGrid.f90
	$(FC) $(FFLAGS) -c type_GaussianGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_Workspace.o: \
type_Workspace.f90
	$(FC) $(FFLAGS) -c type_Workspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_RegularWorkspace.o: \
type_Workspace.o\
type_RegularWorkspace.f90
	$(FC) $(FFLAGS) -c type_RegularWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_GaussianWorkspace.o: \
type_Workspace.o\
type_GaussianWorkspace.f90
	$(FC) $(FFLAGS) -c type_GaussianWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_ThreeDimensionalVector.o: \
type_ThreeDimensionalVector.f90
	$(FC) $(FFLAGS) -c type_ThreeDimensionalVector.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_TrigonometricFunctions.o: \
type_SphericalGrid.o\
type_TrigonometricFunctions.f90
	$(FC) $(FFLAGS) -c type_TrigonometricFunctions.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_SphericalUnitVectors.o: \
type_SphericalGrid.o\
type_TrigonometricFunctions.o\
type_ThreeDimensionalVector.o \
type_SphericalUnitVectors.f90
	$(FC) $(FFLAGS) -c type_SphericalUnitVectors.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_Sphere.o: \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_SphericalGrid.o\
type_Workspace.o\
type_SphericalUnitVectors.o \
type_Sphere.f90
	$(FC) $(FFLAGS) -c type_Sphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_GaussianSphere.o: \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_GaussianGrid.o\
type_GaussianWorkspace.o\
type_SphericalUnitVectors.o \
type_GaussianSphere.f90
	$(FC) $(FFLAGS) -c type_GaussianSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

type_RegularSphere.o: \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_RegularGrid.o\
type_RegularWorkspace.o\
type_SphericalUnitVectors.o \
type_RegularSphere.f90
	$(FC) $(FFLAGS) -c type_RegularSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

spherepack_wrapper_library.o: \
type_ThreeDimensionalVector.o \
type_GaussianSphere.o \
type_RegularSphere.o \
spherepack_wrapper_library.f90
	$(FC) $(FFLAGS) -c spherepack_wrapper_library.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
spherepack_wrapper_library.o \
main.f90
	$(FC) $(FFLAGS) -c main.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

lib: 
	$(AR) rc lib$(LIB_NAME).a $(OBJS)
	mv lib$(LIB_NAME).a $(LIB_DIRECTORY)

all:
	$(MAKE) clean
	$(MAKE) exec; $(EXECUTABLE)
	$(MAKE) lib

clean:
	rm -f $(OBJS) $(LIB) $(LIB_DIRECTORY)/*.mod
	rm -f $(EXECUTABLE)
