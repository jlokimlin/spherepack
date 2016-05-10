
include ../make.inc

LINKS = \
-L$(LIB_DIRECTORY) -lold_spherepack \
-I$(MOD_DIRECTORY)

spherepack_library.o: \
type_ThreeDimensionalVector.o \
type_GaussianSphere.o \
type_RegularSphere.o

type_RegularGrid.o: \
type_SphericalGrid.o

type_GaussianGrid.o: \
type_SphericalGrid.o

type_RegularWorkspace.o: \
type_Workspace.o

type_GaussianWorkspace.o: \
type_Workspace.o

type_TrigonometricFunctions.o: \
type_SphericalGrid.o

type_SphericalUnitVectors.o: \
type_SphericalGrid.o \
type_TrigonometricFunctions.o\
type_ThreeDimensionalVector.o

type_Sphere.o: \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_SphericalGrid.o\
type_Workspace.o\
type_SphericalUnitVectors.o

type_GaussianSphere.o: \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_GaussianGrid.o\
type_GaussianWorkspace.o\
type_SphericalUnitVectors.o

type_RegularSphere.o: \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_RegularGrid.o\
type_RegularWorkspace.o\
type_SphericalUnitVectors.o 

$(LIB): $(OBJS)
	$(AR) -rv $@ $?

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LINKS) -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
clean:
	rm -f $(OBJS) $(LIB) $(MOD_DIRECTORY)/*.mod
