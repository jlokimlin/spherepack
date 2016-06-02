
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(MOD_DIRECTORY)

EXECUTABLE = ./build_lib.exe

geo2math.o: \
geo2math.f90
	$(FC) $(FCFLAGS) -c geo2math.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_FFTpack.o: \
type_FFTpack.f90
	$(FC) $(FCFLAGS) -c type_FFTpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_HFFTpack.o: \
type_HFFTpack.f90
	$(FC) $(FCFLAGS) -c type_HFFTpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SpherepackAux.o: \
type_SpherepackAux.f90
	$(FC) $(FCFLAGS) -c type_SpherepackAux.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 	 
type_LegendreAux.o: \
type_FFTpack.o \
type_LegendreAux.f90
	$(FC) $(FCFLAGS) -c type_LegendreAux.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gaqd.o: \
gaqd.f90
	$(FC) $(FCFLAGS) -c gaqd.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shaec.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shaec.f90
	$(FC) $(FCFLAGS) -c shaec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shaes.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
shaes.f90
	$(FC) $(FCFLAGS) -c shaes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shagc.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
shagc.f90
	$(FC) $(FCFLAGS) -c shagc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shags.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shags.f90
	$(FC) $(FCFLAGS) -c shags.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shsec.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
shsec.f90
	$(FC) $(FCFLAGS) -c shsec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shses.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
shses.f90
	$(FC) $(FCFLAGS) -c shses.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shsgc.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shsgc.f90
	$(FC) $(FCFLAGS) -c shsgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shsgs.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shsgs.f90
	$(FC) $(FCFLAGS) -c shsgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhaec.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
vhaec.f90
	$(FC) $(FCFLAGS) -c vhaec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhaes.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
vhaes.f90
	$(FC) $(FCFLAGS) -c vhaes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhagc.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
vhagc.f90
	$(FC) $(FCFLAGS) -c vhagc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhags.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
vhags.f90
	$(FC) $(FCFLAGS) -c vhags.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsec.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
vhsec.f90
	$(FC) $(FCFLAGS) -c vhsec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhses.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
vhses.f90
	$(FC) $(FCFLAGS) -c vhses.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsgc.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
vhsgc.f90
	$(FC) $(FCFLAGS) -c vhsgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsgs.o: \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
vhsgs.f90
	$(FC) $(FCFLAGS) -c vhsgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divec.o: \
shsec.o \
divec.f90
	$(FC) $(FCFLAGS) -c divec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

dives.o: \
shses.o \
dives.f90
	$(FC) $(FCFLAGS) -c dives.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divgc.o: \
shsgc.o \
divgc.f90
	$(FC) $(FCFLAGS) -c divgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divgs.o: \
shsgs.o \
divgs.f90
	$(FC) $(FCFLAGS) -c divgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SphericalGrid.o: \
type_SphericalGrid.f90
	$(FC) $(FCFLAGS) -c type_SphericalGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
type_RegularGrid.o: \
type_SphericalGrid.o\
type_RegularGrid.f90
	$(FC) $(FCFLAGS) -c type_RegularGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
type_GaussianGrid.o: \
type_SphericalGrid.o \
gaqd.o \
type_GaussianGrid.f90
	$(FC) $(FCFLAGS) -c type_GaussianGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Workspace.o: \
type_Workspace.f90
	$(FC) $(FCFLAGS) -c type_Workspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RegularWorkspace.o: \
type_Workspace.o\
shaes.o \
shses.o \
vhaes.o \
vhses.o \
type_RegularWorkspace.f90
	$(FC) $(FCFLAGS) -c type_RegularWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_GaussianWorkspace.o: \
type_Workspace.o\
shags.o \
shsgs.o \
vhags.o \
vhsgs.o \
type_GaussianWorkspace.f90
	$(FC) $(FCFLAGS) -c type_GaussianWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_ThreeDimensionalVector.o: \
type_ThreeDimensionalVector.f90
	$(FC) $(FCFLAGS) -c type_ThreeDimensionalVector.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_TrigonometricFunctions.o: \
type_SphericalGrid.o\
type_TrigonometricFunctions.f90
	$(FC) $(FCFLAGS) -c type_TrigonometricFunctions.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SphericalUnitVectors.o: \
type_SphericalGrid.o\
type_TrigonometricFunctions.o\
type_ThreeDimensionalVector.o \
type_SphericalUnitVectors.f90
	$(FC) $(FCFLAGS) -c type_SphericalUnitVectors.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Sphere.o: \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_SphericalGrid.o\
type_Workspace.o\
type_SphericalUnitVectors.o \
type_Sphere.f90
	$(FC) $(FCFLAGS) -c type_Sphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_GaussianSphere.o: \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_GaussianGrid.o\
type_GaussianWorkspace.o\
type_SphericalUnitVectors.o \
shags.o \
shsgs.o \
vhags.o \
vhsgs.o \
type_GaussianSphere.f90
	$(FC) $(FCFLAGS) -c type_GaussianSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RegularSphere.o: \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_RegularGrid.o\
type_RegularWorkspace.o\
type_SphericalUnitVectors.o \
shaes.o \
shses.o \
vhaes.o \
vhses.o \
type_RegularSphere.f90
	$(FC) $(FCFLAGS) -c type_RegularSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

spherepack_library.o: \
type_ThreeDimensionalVector.o \
type_GaussianSphere.o \
type_RegularSphere.o \
spherepack_library.f90
	$(FC) $(FCFLAGS) -c spherepack_library.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

main.o: \
spherepack_library.o \
main.f90
	$(FC) $(FCFLAGS) -c main.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FCFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

lib: 
	$(AR) rc lib$(LIB_NAME).a $(OBJ_DIRECTORY)/*.o
	mv lib$(LIB_NAME).a $(LIB_DIRECTORY)

all:
	$(MAKE) clean
	$(MAKE) exec; $(EXECUTABLE)
	$(MAKE) lib

clean:
	rm -f $(OBJS) $(LIB) $(MOD_DIRECTORY)/*.mod
	rm -f $(EXECUTABLE)