
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(MOD_DIRECTORY)

EXECUTABLE = ./build_lib.exe

spherepack_precision.o: \
spherepack_precision.f90
	$(FC) $(FFLAGS) -c spherepack_precision.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

geo2math.o: \
spherepack_precision.o \
geo2math.f90
	$(FC) $(FFLAGS) -c geo2math.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ihgeod.o: \
spherepack_precision.o \
ihgeod.f90
	$(FC) $(FFLAGS) -c ihgeod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_FFTpack.o: \
spherepack_precision.o \
type_FFTpack.f90
	$(FC) $(FFLAGS) -c type_FFTpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_HFFTpack.o: \
spherepack_precision.o \
type_HFFTpack.f90
	$(FC) $(FFLAGS) -c type_HFFTpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SpherepackAux.o: \
spherepack_precision.o \
type_SpherepackAux.f90
	$(FC) $(FFLAGS) -c type_SpherepackAux.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 	 
type_LegendreAux.o: \
spherepack_precision.o \
type_FFTpack.o \
type_LegendreAux.f90
	$(FC) $(FFLAGS) -c type_LegendreAux.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shpg.o: \
spherepack_precision.o \
type_HFFTpack.o \
shpg.f90
	$(FC) $(FFLAGS) -c shpg.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

shpe.o: \
spherepack_precision.o \
type_HFFTpack.o \
shpe.f90
	$(FC) $(FFLAGS) -c shpe.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gaussian_latitudes_and_weights_routines.o: \
spherepack_precision.o \
gaussian_latitudes_and_weights_routines.f90
	$(FC) $(FFLAGS) -c gaussian_latitudes_and_weights_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_analysis_routines.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaussian_latitudes_and_weights_routines.o\
scalar_analysis_routines.f90
	$(FC) $(FFLAGS) -c scalar_analysis_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_analysis_regular_grid.o: \
scalar_analysis_routines.o \
scalar_analysis_regular_grid.f90
	$(FC) $(FFLAGS) -c scalar_analysis_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_analysis_regular_grid_saved.o: \
scalar_analysis_routines.o \
scalar_analysis_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c scalar_analysis_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_analysis_gaussian_grid.o: \
scalar_analysis_routines.o \
scalar_analysis_gaussian_grid.f90
	$(FC) $(FFLAGS) -c scalar_analysis_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_analysis_gaussian_grid_saved.o: \
scalar_analysis_routines.o \
scalar_analysis_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c scalar_analysis_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_synthesis_routines.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaussian_latitudes_and_weights_routines.o\
scalar_synthesis_routines.f90
	$(FC) $(FFLAGS) -c scalar_synthesis_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_synthesis_regular_grid.o: \
scalar_synthesis_routines.o \
scalar_analysis_regular_grid.f90
	$(FC) $(FFLAGS) -c scalar_synthesis_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_synthesis_regular_grid_saved.o: \
scalar_synthesis_routines.o \
scalar_synthesis_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c scalar_synthesis_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_synthesis_gaussian_grid.o: \
scalar_synthesis_routines.o \
scalar_synthesis_gaussian_grid.f90
	$(FC) $(FFLAGS) -c scalar_synthesis_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_synthesis_gaussian_grid_saved.o: \
scalar_synthesis_routines.o \
scalar_synthesis_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c scalar_synthesis_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_analysis_routines.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaussian_latitudes_and_weights_routines.o\
vector_analysis_routines.f90
	$(FC) $(FFLAGS) -c vector_analysis_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_analysis_regular_grid.o: \
vector_analysis_routines.o \
vector_analysis_regular_grid.f90
	$(FC) $(FFLAGS) -c vector_analysis_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_analysis_regular_grid_saved.o: \
vector_analysis_routines.o \
vector_analysis_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c vector_analysis_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_analysis_gaussian_grid.o: \
vector_analysis_routines.o \
vector_analysis_gaussian_grid.f90
	$(FC) $(FFLAGS) -c vector_analysis_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_analysis_gaussian_grid_saved.o: \
vector_analysis_routines.o \
vector_analysis_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c vector_analysis_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsec.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vhsec.f90
	$(FC) $(FFLAGS) -c vhsec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhses.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vhses.f90
	$(FC) $(FFLAGS) -c vhses.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsgc.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaussian_latitudes_and_weights_routines.o \
vhsgc.f90
	$(FC) $(FFLAGS) -c vhsgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vhsgs.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaussian_latitudes_and_weights_routines.o \
vhsgs.f90
	$(FC) $(FFLAGS) -c vhsgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtsec.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vtsec.f90
	$(FC) $(FFLAGS) -c vtsec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtses.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
vtses.f90
	$(FC) $(FFLAGS) -c vtses.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtsgc.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaussian_latitudes_and_weights_routines.o \
vtsgc.f90
	$(FC) $(FFLAGS) -c vtsgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vtsgs.o: \
spherepack_precision.o \
type_HFFTpack.o \
type_SpherepackAux.o\
gaussian_latitudes_and_weights_routines.o \
vtsgs.f90
	$(FC) $(FFLAGS) -c vtsgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)


sshifte.o: \
spherepack_precision.o \
type_HFFTpack.o \
sshifte.f90
	$(FC) $(FFLAGS) -c sshifte.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vshifte.o: \
spherepack_precision.o \
type_HFFTpack.o \
vshifte.f90
	$(FC) $(FFLAGS) -c vshifte.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

trssph.o: \
spherepack_precision.o \
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
trssph.f90
	$(FC) $(FFLAGS) -c trssph.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

trvsph.o: \
spherepack_precision.o \
type_HFFTpack.o \
vector_analysis_routines.o \
vhsec.o \
vector_analysis_routines.o \
vhsgc.o \
trvsph.f90
	$(FC) $(FFLAGS) -c trvsph.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradec.o: \
spherepack_precision.o \
vhsec.o \
gradec.f90
	$(FC) $(FFLAGS) -c gradec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

grades.o: \
spherepack_precision.o \
vhses.o \
grades.f90
	$(FC) $(FFLAGS) -c grades.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradgc.o: \
spherepack_precision.o \
vhsgc.o \
gradgc.f90
	$(FC) $(FFLAGS) -c gradgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradgs.o: \
spherepack_precision.o \
vhsgs.o \
gradgs.f90
	$(FC) $(FFLAGS) -c gradgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igradec.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
igradec.f90
	$(FC) $(FFLAGS) -c igradec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igrades.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
igrades.f90
	$(FC) $(FFLAGS) -c igrades.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igradgc.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
igradgc.f90
	$(FC) $(FFLAGS) -c igradgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

igradgs.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
igradgs.f90
	$(FC) $(FFLAGS) -c igradgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divergence_routines.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
divergence_routines.f90
	$(FC) $(FFLAGS) -c divergence_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divergence_regular_grid.o: \
divergence_routines.o \
divergence_regular_grid.f90
	$(FC) $(FFLAGS) -c divergence_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divergence_regular_grid_saved.o: \
divergence_routines.o \
divergence_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c divergence_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divergence_gaussian_grid.o: \
divergence_routines.o \
divergence_gaussian_grid.f90
	$(FC) $(FFLAGS) -c divergence_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divergence_gaussian_grid_saved.o: \
divergence_routines.o \
divergence_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c divergence_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idivec.o: \
spherepack_precision.o \
vhsec.o \
idivec.f90
	$(FC) $(FFLAGS) -c idivec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idives.o: \
spherepack_precision.o \
vhses.o \
idives.f90
	$(FC) $(FFLAGS) -c idives.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idivgc.o: \
spherepack_precision.o \
vhsgc.o \
idivgc.f90
	$(FC) $(FFLAGS) -c idivgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idivgs.o: \
spherepack_precision.o \
vhsgs.o \
idivgs.f90
	$(FC) $(FFLAGS) -c idivgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtec.o: \
spherepack_precision.o \
vhsec.o \
idvtec.f90
	$(FC) $(FFLAGS) -c idvtec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtes.o: \
spherepack_precision.o \
vhses.o \
idvtes.f90
	$(FC) $(FFLAGS) -c idvtes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtgc.o: \
spherepack_precision.o \
vhsgc.o \
idvtgc.f90
	$(FC) $(FFLAGS) -c idvtgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtgs.o: \
spherepack_precision.o \
vhsgs.o \
idvtgs.f90
	$(FC) $(FFLAGS) -c idvtgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapec.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
slapec.f90
	$(FC) $(FFLAGS) -c slapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapes.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
slapes.f90
	$(FC) $(FFLAGS) -c slapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapgc.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
slapgc.f90
	$(FC) $(FFLAGS) -c slapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

slapgs.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
slapgs.f90
	$(FC) $(FFLAGS) -c slapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapec.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
islapec.f90
	$(FC) $(FFLAGS) -c islapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapes.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
islapes.f90
	$(FC) $(FFLAGS) -c islapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapgc.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
islapgc.f90
	$(FC) $(FFLAGS) -c islapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

islapgs.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
islapgs.f90
	$(FC) $(FFLAGS) -c islapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpec.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
sfvpec.f90
	$(FC) $(FFLAGS) -c sfvpec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpes.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
sfvpes.f90
	$(FC) $(FFLAGS) -c sfvpes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpgc.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
sfvpgc.f90
	$(FC) $(FFLAGS) -c sfvpgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sfvpgs.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
sfvpgs.f90
	$(FC) $(FFLAGS) -c sfvpgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpec.o: \
spherepack_precision.o \
vhsec.o \
isfvpec.f90
	$(FC) $(FFLAGS) -c isfvpec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpes.o: \
spherepack_precision.o \
vhses.o \
isfvpes.f90
	$(FC) $(FFLAGS) -c isfvpes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpgc.o: \
spherepack_precision.o \
vhsgc.o \
isfvpgc.f90
	$(FC) $(FFLAGS) -c isfvpgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpgs.o: \
spherepack_precision.o \
vhsgs.o \
isfvpgs.f90
	$(FC) $(FFLAGS) -c isfvpgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtec.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
vrtec.f90
	$(FC) $(FFLAGS) -c vrtec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtes.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
vrtes.f90
	$(FC) $(FFLAGS) -c vrtes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtgc.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
vrtgc.f90
	$(FC) $(FFLAGS) -c vrtgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vrtgs.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
vrtgs.f90
	$(FC) $(FFLAGS) -c vrtgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtec.o: \
spherepack_precision.o \
vhsec.o \
ivrtec.f90
	$(FC) $(FFLAGS) -c ivrtec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtes.o: \
spherepack_precision.o \
vhses.o \
ivrtes.f90
	$(FC) $(FFLAGS) -c ivrtes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtgc.o: \
spherepack_precision.o \
vhsgc.o \
ivrtgc.f90
	$(FC) $(FFLAGS) -c ivrtgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivrtgs.o: \
spherepack_precision.o \
vhsgs.o \
ivrtgs.f90
	$(FC) $(FFLAGS) -c ivrtgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapec.o: \
spherepack_precision.o \
vhsec.o \
vlapec.f90
	$(FC) $(FFLAGS) -c vlapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapes.o: \
spherepack_precision.o \
vhses.o \
vlapes.f90
	$(FC) $(FFLAGS) -c vlapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapgc.o: \
spherepack_precision.o \
vhsgc.o \
vlapgc.f90
	$(FC) $(FFLAGS) -c vlapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vlapgs.o: \
spherepack_precision.o \
vhsgs.o \
vlapgs.f90
	$(FC) $(FFLAGS) -c vlapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapec.o: \
spherepack_precision.o \
vhsec.o \
ivlapec.f90
	$(FC) $(FFLAGS) -c ivlapec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapes.o: \
spherepack_precision.o \
vhses.o \
ivlapes.f90
	$(FC) $(FFLAGS) -c ivlapes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapgc.o: \
spherepack_precision.o \
vhsgc.o \
ivlapgc.f90
	$(FC) $(FFLAGS) -c ivlapgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

ivlapgs.o: \
spherepack_precision.o \
vhsgs.o \
ivlapgs.f90
	$(FC) $(FFLAGS) -c ivlapgs.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SphericalGrid.o: \
spherepack_precision.o \
type_SphericalGrid.f90
	$(FC) $(FFLAGS) -c type_SphericalGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
type_RegularGrid.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_RegularGrid.f90
	$(FC) $(FFLAGS) -c type_RegularGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 
type_GaussianGrid.o: \
spherepack_precision.o \
type_SphericalGrid.o \
gaussian_latitudes_and_weights_routines.o \
type_GaussianGrid.f90
	$(FC) $(FFLAGS) -c type_GaussianGrid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Workspace.o: \
spherepack_precision.o \
type_Workspace.f90
	$(FC) $(FFLAGS) -c type_Workspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RegularWorkspace.o: \
spherepack_precision.o \
type_Workspace.o\
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vhses.o \
type_RegularWorkspace.f90
	$(FC) $(FFLAGS) -c type_RegularWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_GaussianWorkspace.o: \
spherepack_precision.o \
type_Workspace.o\
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vhsgs.o \
type_GaussianWorkspace.f90
	$(FC) $(FFLAGS) -c type_GaussianWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Vector3D.o: \
spherepack_precision.o \
type_Vector3D.f90
	$(FC) $(FFLAGS) -c type_Vector3D.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_TrigonometricFunctions.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_TrigonometricFunctions.f90
	$(FC) $(FFLAGS) -c type_TrigonometricFunctions.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SphericalUnitVectors.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_TrigonometricFunctions.o\
type_Vector3D.o \
type_SphericalUnitVectors.f90
	$(FC) $(FFLAGS) -c type_SphericalUnitVectors.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Sphere.o: \
spherepack_precision.o \
type_TrigonometricFunctions.o \
type_Vector3D.o \
type_SphericalGrid.o\
type_Workspace.o\
type_SphericalUnitVectors.o \
type_Sphere.f90
	$(FC) $(FFLAGS) -c type_Sphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_GaussianSphere.o: \
spherepack_precision.o \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_Vector3D.o \
type_GaussianGrid.o\
type_GaussianWorkspace.o\
type_SphericalUnitVectors.o \
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vhsgs.o \
type_GaussianSphere.f90
	$(FC) $(FFLAGS) -c type_GaussianSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RegularSphere.o: \
spherepack_precision.o \
type_Sphere.o \
type_TrigonometricFunctions.o \
type_Vector3D.o \
type_RegularGrid.o\
type_RegularWorkspace.o\
type_SphericalUnitVectors.o \
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vhses.o \
type_RegularSphere.f90
	$(FC) $(FFLAGS) -c type_RegularSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

spherepack_library.o: \
spherepack_precision.o \
divergence_routines.o \
gaussian_latitudes_and_weights_routines.o \
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
scalar_analysis_routines.o \
shpe.o \
shpg.o \
scalar_synthesis_routines.o \
slapec.o \
slapes.o \
slapgc.o \
slapgs.o \
sshifte.o \
trssph.o \
trvsph.o \
vector_analysis_routines.o \
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
type_Vector3D.o \
type_GaussianSphere.o \
type_RegularSphere.o \
spherepack_library.f90
	$(FC) $(FFLAGS) -c spherepack_library.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

main.o: \
spherepack_library.o \
main.f90
	$(FC) $(FFLAGS) -c main.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

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