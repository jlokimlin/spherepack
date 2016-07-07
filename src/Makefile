
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(MOD_DIRECTORY)

EXECUTABLE = ./build_lib.exe

spherepack_precision.o: \
spherepack_precision.f90
	$(FC) $(FCFLAGS) -c spherepack_precision.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

geo2math.o: \
spherepack_precision.o \
geo2math.f90
	$(FC) $(FCFLAGS) -c geo2math.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ihgeod.o: \
spherepack_precision.o \
ihgeod.f90
	$(FC) $(FCFLAGS) -c ihgeod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_FFTpack.o: \
spherepack_precision.o \
type_FFTpack.f90
	$(FC) $(FCFLAGS) -c type_FFTpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_HFFTpack.o: \
spherepack_precision.o \
type_HFFTpack.f90
	$(FC) $(FCFLAGS) -c type_HFFTpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SpherepackAux.o: \
spherepack_precision.o \
type_SpherepackAux.f90
	$(FC) $(FCFLAGS) -c type_SpherepackAux.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 	 
type_LegendreAux.o: \
spherepack_precision.o \
type_FFTpack.o \
type_LegendreAux.f90
	$(FC) $(FCFLAGS) -c type_LegendreAux.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shpg.o: \
spherepack_precision.o \
type_HFFTpack.o \
shpg.f90
	$(FC) $(FCFLAGS) -c shpg.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shpe.o: \
spherepack_precision.o \
type_HFFTpack.o \
shpe.f90
	$(FC) $(FCFLAGS) -c shpe.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gaqd.o: \
spherepack_precision.o \
gaqd.f90
	$(FC) $(FCFLAGS) -c gaqd.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shaec.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shaec.f90
	$(FC) $(FCFLAGS) -c shaec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shaes.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
shaes.f90
	$(FC) $(FCFLAGS) -c shaes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shagc.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
shagc.f90
	$(FC) $(FCFLAGS) -c shagc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shags.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shags.f90
	$(FC) $(FCFLAGS) -c shags.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shsec.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
shsec.f90
	$(FC) $(FCFLAGS) -c shsec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shses.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
shses.f90
	$(FC) $(FCFLAGS) -c shses.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shsgc.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shsgc.f90
	$(FC) $(FCFLAGS) -c shsgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shsgs.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
shsgs.f90
	$(FC) $(FCFLAGS) -c shsgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhaec.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vhaec.f90
	$(FC) $(FCFLAGS) -c vhaec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhaes.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vhaes.f90
	$(FC) $(FCFLAGS) -c vhaes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhagc.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
vhagc.f90
	$(FC) $(FCFLAGS) -c vhagc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhags.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o\
vhags.f90
	$(FC) $(FCFLAGS) -c vhags.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsec.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vhsec.f90
	$(FC) $(FCFLAGS) -c vhsec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhses.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vhses.f90
	$(FC) $(FCFLAGS) -c vhses.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsgc.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
vhsgc.f90
	$(FC) $(FCFLAGS) -c vhsgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsgs.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
vhsgs.f90
	$(FC) $(FCFLAGS) -c vhsgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtsec.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vtsec.f90
	$(FC) $(FCFLAGS) -c vtsec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtses.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vtses.f90
	$(FC) $(FCFLAGS) -c vtses.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtsgc.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
vtsgc.f90
	$(FC) $(FCFLAGS) -c vtsgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtsgs.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaqd.o \
vtsgs.f90
	$(FC) $(FCFLAGS) -c vtsgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)


sshifte.o: \
spherepack_precision.o \
type_HFFTpack.o \
sshifte.f90
	$(FC) $(FCFLAGS) -c sshifte.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vshifte.o: \
spherepack_precision.o \
type_HFFTpack.o \
vshifte.f90
	$(FC) $(FCFLAGS) -c vshifte.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

trssph.o: \
spherepack_precision.o \
type_HFFTpack.o \
shaec.o \
shsec.o \
shagc.o \
shsgc.o \
trssph.f90
	$(FC) $(FCFLAGS) -c trssph.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

trvsph.o: \
spherepack_precision.o \
type_HFFTpack.o \
vhaec.o \
vhsec.o \
vhagc.o \
vhsgc.o \
trvsph.f90
	$(FC) $(FCFLAGS) -c trvsph.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradec.o: \
spherepack_precision.o \
vhsec.o \
gradec.f90
	$(FC) $(FCFLAGS) -c gradec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

grades.o: \
spherepack_precision.o \
vhses.o \
grades.f90
	$(FC) $(FCFLAGS) -c grades.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradgc.o: \
spherepack_precision.o \
vhsgc.o \
gradgc.f90
	$(FC) $(FCFLAGS) -c gradgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradgs.o: \
spherepack_precision.o \
vhsgs.o \
gradgs.f90
	$(FC) $(FCFLAGS) -c gradgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igradec.o: \
spherepack_precision.o \
shsec.o \
igradec.f90
	$(FC) $(FCFLAGS) -c igradec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igrades.o: \
spherepack_precision.o \
shses.o \
igrades.f90
	$(FC) $(FCFLAGS) -c igrades.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igradgc.o: \
spherepack_precision.o \
shsgc.o \
igradgc.f90
	$(FC) $(FCFLAGS) -c igradgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igradgs.o: \
spherepack_precision.o \
shsgs.o \
igradgs.f90
	$(FC) $(FCFLAGS) -c igradgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divec.o: \
spherepack_precision.o \
shsec.o \
divec.f90
	$(FC) $(FCFLAGS) -c divec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

dives.o: \
spherepack_precision.o \
shses.o \
dives.f90
	$(FC) $(FCFLAGS) -c dives.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divgc.o: \
spherepack_precision.o \
shsgc.o \
divgc.f90
	$(FC) $(FCFLAGS) -c divgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divgs.o: \
spherepack_precision.o \
shsgs.o \
divgs.f90
	$(FC) $(FCFLAGS) -c divgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idivec.o: \
spherepack_precision.o \
vhsec.o \
idivec.f90
	$(FC) $(FCFLAGS) -c idivec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idives.o: \
spherepack_precision.o \
vhses.o \
idives.f90
	$(FC) $(FCFLAGS) -c idives.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idivgc.o: \
spherepack_precision.o \
vhsgc.o \
idivgc.f90
	$(FC) $(FCFLAGS) -c idivgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idivgs.o: \
spherepack_precision.o \
vhsgs.o \
idivgs.f90
	$(FC) $(FCFLAGS) -c idivgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtec.o: \
spherepack_precision.o \
vhsec.o \
idvtec.f90
	$(FC) $(FCFLAGS) -c idvtec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtes.o: \
spherepack_precision.o \
vhses.o \
idvtes.f90
	$(FC) $(FCFLAGS) -c idvtes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtgc.o: \
spherepack_precision.o \
vhsgc.o \
idvtgc.f90
	$(FC) $(FCFLAGS) -c idvtgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtgs.o: \
spherepack_precision.o \
vhsgs.o \
idvtgs.f90
	$(FC) $(FCFLAGS) -c idvtgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapec.o: \
spherepack_precision.o \
shsec.o \
slapec.f90
	$(FC) $(FCFLAGS) -c slapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapes.o: \
spherepack_precision.o \
shses.o \
slapes.f90
	$(FC) $(FCFLAGS) -c slapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapgc.o: \
spherepack_precision.o \
shsgc.o \
slapgc.f90
	$(FC) $(FCFLAGS) -c slapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapgs.o: \
spherepack_precision.o \
shsgs.o \
slapgs.f90
	$(FC) $(FCFLAGS) -c slapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapec.o: \
spherepack_precision.o \
shsec.o \
islapec.f90
	$(FC) $(FCFLAGS) -c islapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapes.o: \
spherepack_precision.o \
shses.o \
islapes.f90
	$(FC) $(FCFLAGS) -c islapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapgc.o: \
spherepack_precision.o \
shsgc.o \
islapgc.f90
	$(FC) $(FCFLAGS) -c islapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapgs.o: \
spherepack_precision.o \
shsgs.o \
islapgs.f90
	$(FC) $(FCFLAGS) -c islapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpec.o: \
spherepack_precision.o \
shsec.o \
sfvpec.f90
	$(FC) $(FCFLAGS) -c sfvpec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpes.o: \
spherepack_precision.o \
shses.o \
sfvpes.f90
	$(FC) $(FCFLAGS) -c sfvpes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpgc.o: \
spherepack_precision.o \
shsgc.o \
sfvpgc.f90
	$(FC) $(FCFLAGS) -c sfvpgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpgs.o: \
spherepack_precision.o \
shsgs.o \
sfvpgs.f90
	$(FC) $(FCFLAGS) -c sfvpgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpec.o: \
spherepack_precision.o \
vhsec.o \
isfvpec.f90
	$(FC) $(FCFLAGS) -c isfvpec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpes.o: \
spherepack_precision.o \
vhses.o \
isfvpes.f90
	$(FC) $(FCFLAGS) -c isfvpes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpgc.o: \
spherepack_precision.o \
vhsgc.o \
isfvpgc.f90
	$(FC) $(FCFLAGS) -c isfvpgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpgs.o: \
spherepack_precision.o \
vhsgs.o \
isfvpgs.f90
	$(FC) $(FCFLAGS) -c isfvpgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtec.o: \
spherepack_precision.o \
shsec.o \
vrtec.f90
	$(FC) $(FCFLAGS) -c vrtec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtes.o: \
spherepack_precision.o \
shses.o \
vrtes.f90
	$(FC) $(FCFLAGS) -c vrtes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtgc.o: \
spherepack_precision.o \
shsgc.o \
vrtgc.f90
	$(FC) $(FCFLAGS) -c vrtgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtgs.o: \
spherepack_precision.o \
shsgs.o \
vrtgs.f90
	$(FC) $(FCFLAGS) -c vrtgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtec.o: \
spherepack_precision.o \
vhsec.o \
ivrtec.f90
	$(FC) $(FCFLAGS) -c ivrtec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtes.o: \
spherepack_precision.o \
vhses.o \
ivrtes.f90
	$(FC) $(FCFLAGS) -c ivrtes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtgc.o: \
spherepack_precision.o \
vhsgc.o \
ivrtgc.f90
	$(FC) $(FCFLAGS) -c ivrtgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtgs.o: \
spherepack_precision.o \
vhsgs.o \
ivrtgs.f90
	$(FC) $(FCFLAGS) -c ivrtgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapec.o: \
spherepack_precision.o \
vhsec.o \
vlapec.f90
	$(FC) $(FCFLAGS) -c vlapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapes.o: \
spherepack_precision.o \
vhses.o \
vlapes.f90
	$(FC) $(FCFLAGS) -c vlapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapgc.o: \
spherepack_precision.o \
vhsgc.o \
vlapgc.f90
	$(FC) $(FCFLAGS) -c vlapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapgs.o: \
spherepack_precision.o \
vhsgs.o \
vlapgs.f90
	$(FC) $(FCFLAGS) -c vlapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapec.o: \
spherepack_precision.o \
vhsec.o \
ivlapec.f90
	$(FC) $(FCFLAGS) -c ivlapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapes.o: \
spherepack_precision.o \
vhses.o \
ivlapes.f90
	$(FC) $(FCFLAGS) -c ivlapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapgc.o: \
spherepack_precision.o \
vhsgc.o \
ivlapgc.f90
	$(FC) $(FCFLAGS) -c ivlapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapgs.o: \
spherepack_precision.o \
vhsgs.o \
ivlapgs.f90
	$(FC) $(FCFLAGS) -c ivlapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SphericalGrid.o: \
spherepack_precision.o \
type_SphericalGrid.f90
	$(FC) $(FCFLAGS) -c type_SphericalGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
type_RegularGrid.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_RegularGrid.f90
	$(FC) $(FCFLAGS) -c type_RegularGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
type_GaussianGrid.o: \
spherepack_precision.o \
type_SphericalGrid.o \
gaqd.o \
type_GaussianGrid.f90
	$(FC) $(FCFLAGS) -c type_GaussianGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Workspace.o: \
spherepack_precision.o \
type_Workspace.f90
	$(FC) $(FCFLAGS) -c type_Workspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RegularWorkspace.o: \
spherepack_precision.o \
type_Workspace.o\
shaes.o \
shses.o \
vhaes.o \
vhses.o \
type_RegularWorkspace.f90
	$(FC) $(FCFLAGS) -c type_RegularWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_GaussianWorkspace.o: \
spherepack_precision.o \
type_Workspace.o\
shags.o \
shsgs.o \
vhags.o \
vhsgs.o \
type_GaussianWorkspace.f90
	$(FC) $(FCFLAGS) -c type_GaussianWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_ThreeDimensionalVector.o: \
spherepack_precision.o \
type_ThreeDimensionalVector.f90
	$(FC) $(FCFLAGS) -c type_ThreeDimensionalVector.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_TrigonometricFunctions.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_TrigonometricFunctions.f90
	$(FC) $(FCFLAGS) -c type_TrigonometricFunctions.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SphericalUnitVectors.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_TrigonometricFunctions.o\
type_ThreeDimensionalVector.o \
type_SphericalUnitVectors.f90
	$(FC) $(FCFLAGS) -c type_SphericalUnitVectors.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Sphere.o: \
spherepack_precision.o \
type_TrigonometricFunctions.o \
type_ThreeDimensionalVector.o \
type_SphericalGrid.o\
type_Workspace.o\
type_SphericalUnitVectors.o \
type_Sphere.f90
	$(FC) $(FCFLAGS) -c type_Sphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_GaussianSphere.o: \
spherepack_precision.o \
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
spherepack_precision.o \
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
spherepack_precision.o \
divec.o \
dives.o \
divgc.o \
divgs.o \
gaqd.o \
geo2math.o \
gradec.o \
grades.o \
gradgc.o \
gradgs.o \
idivec.o \
idives.o \
idivgc.o \
idivgs.o \
idvtec.o \
idvtes.o \
idvtgc.o \
idvtgs.o \
igradec.o \
igrades.o \
igradgc.o \
igradgs.o \
ihgeod.o \
isfvpec.o \
isfvpes.o \
isfvpgc.o \
isfvpgs.o \
islapec.o \
islapes.o \
islapgc.o \
islapgs.o \
ivlapec.o \
ivlapes.o \
ivlapgc.o \
ivlapgs.o \
ivrtec.o \
ivrtes.o \
ivrtgc.o \
ivrtgs.o \
sfvpec.o \
sfvpes.o \
sfvpgc.o \
sfvpgs.o \
shaec.o \
shaes.o \
shagc.o \
shags.o \
shpe.o \
shpg.o \
shsec.o \
shses.o \
shsgc.o \
shsgs.o \
slapec.o \
slapes.o \
slapgc.o \
slapgs.o \
sshifte.o \
trssph.o \
trvsph.o \
vhaec.o \
vhaes.o \
vhagc.o \
vhags.o \
vhsec.o \
vhses.o \
vhsgc.o \
vhsgs.o \
vlapec.o \
vlapes.o \
vlapgc.o \
vlapgs.o \
vrtec.o \
vrtes.o \
vrtgc.o \
vrtgs.o \
vshifte.o \
vtsec.o \
vtses.o \
vtsgc.o \
vtsgs.o \
type_HFFTpack.o \
type_FFTpack.o \
type_LegendreAux.o \
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