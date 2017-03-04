
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(MOD_DIRECTORY)

EXECUTABLE = ./build_lib.exe

spherepack_precision.o: \
spherepack_precision.f90
	$(FC) $(FFLAGS) -c spherepack_precision.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

coordinate_transfer_routines.o: \
spherepack_precision.o \
coordinate_transfer_routines.f90
	$(FC) $(FFLAGS) -c coordinate_transfer_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

icosahedral_geodesic_routines.o: \
spherepack_precision.o \
icosahedral_geodesic_routines.f90
	$(FC) $(FFLAGS) -c icosahedral_geodesic_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_FFTpack.o: \
spherepack_precision.o \
type_FFTpack.f90
	$(FC) $(FFLAGS) -c type_FFTpack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RealPeriodicFastFourierTransform.o: \
spherepack_precision.o \
type_RealPeriodicFastFourierTransform.f90
	$(FC) $(FFLAGS) -c type_RealPeriodicFastFourierTransform.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SpherepackUtility.o: \
type_RealPeriodicFastFourierTransform.o \
type_SpherepackUtility.f90
	$(FC) $(FFLAGS) -c type_SpherepackUtility.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)
	 	 
type_LegendreAux.o: \
spherepack_precision.o \
type_FFTpack.o \
type_LegendreAux.f90
	$(FC) $(FFLAGS) -c type_LegendreAux.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)


gaussian_latitudes_and_weights_routines.o: \
spherepack_precision.o \
gaussian_latitudes_and_weights_routines.f90
	$(FC) $(FFLAGS) -c gaussian_latitudes_and_weights_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_projection_routines.o: \
spherepack_precision.o \
type_SpherepackUtility.o \
gaussian_latitudes_and_weights_routines.o \
scalar_projection_routines.f90
	$(FC) $(FFLAGS) -c scalar_projection_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_projection_regular_grid.o: \
scalar_projection_routines.o \
scalar_projection_regular_grid.f90
	$(FC) $(FFLAGS) -c scalar_projection_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_projection_gaussian_grid.o: \
scalar_projection_routines.o \
scalar_projection_gaussian_grid.f90
	$(FC) $(FFLAGS) -c scalar_projection_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_analysis_routines.o: \
spherepack_precision.o \
type_SpherepackUtility.o\
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
type_SpherepackUtility.o\
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
type_SpherepackUtility.o\
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

vector_synthesis_routines.o: \
spherepack_precision.o \
type_SpherepackUtility.o\
gaussian_latitudes_and_weights_routines.o \
vector_synthesis_routines.f90
	$(FC) $(FFLAGS) -c vector_synthesis_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_synthesis_regular_grid.o: \
vector_synthesis_routines.o \
vector_synthesis_regular_grid.f90
	$(FC) $(FFLAGS) -c vector_synthesis_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_synthesis_regular_grid_saved.o: \
vector_synthesis_routines.o \
vector_synthesis_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c vector_synthesis_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_synthesis_gaussian_grid.o: \
vector_synthesis_routines.o \
vector_synthesis_gaussian_grid.f90
	$(FC) $(FFLAGS) -c vector_synthesis_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_synthesis_gaussian_grid_saved.o: \
vector_synthesis_routines.o \
vector_synthesis_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c vector_synthesis_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

colatitudinal_derivative_routines.o: \
spherepack_precision.o \
type_SpherepackUtility.o\
gaussian_latitudes_and_weights_routines.o \
colatitudinal_derivative_routines.f90
	$(FC) $(FFLAGS) -c colatitudinal_derivative_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

colatitudinal_derivative_regular_grid.o: \
colatitudinal_derivative_routines.o \
colatitudinal_derivative_regular_grid.f90
	$(FC) $(FFLAGS) -c colatitudinal_derivative_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

colatitudinal_derivative_regular_grid_saved.o: \
colatitudinal_derivative_routines.o \
colatitudinal_derivative_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c colatitudinal_derivative_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

colatitudinal_derivative_gaussian_grid.o: \
colatitudinal_derivative_routines.o \
colatitudinal_derivative_gaussian_grid.f90
	$(FC) $(FFLAGS) -c colatitudinal_derivative_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

colatitudinal_derivative_gaussian_grid_saved.o: \
colatitudinal_derivative_routines.o \
colatitudinal_derivative_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c colatitudinal_derivative_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

sshifte.o: \
spherepack_precision.o \
type_RealPeriodicFastFourierTransform.o \
sshifte.f90
	$(FC) $(FFLAGS) -c sshifte.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vshifte.o: \
spherepack_precision.o \
type_RealPeriodicFastFourierTransform.o \
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
type_RealPeriodicFastFourierTransform.o \
vector_analysis_routines.o \
vector_synthesis_routines.o \
vector_analysis_routines.o \
vector_synthesis_routines.o \
trvsph.f90
	$(FC) $(FFLAGS) -c trvsph.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradient_routines.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
vector_synthesis_routines.o \
gradient_routines.f90
	$(FC) $(FFLAGS) -c gradient_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradient_regular_grid.o: \
gradient_routines.o \
gradient_regular_grid.f90
	$(FC) $(FFLAGS) -c gradient_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradient_regular_grid_saved.o: \
gradient_routines.o \
gradient_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c gradient_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradient_gaussian_grid.o: \
gradient_routines.o \
gradient_gaussian_grid.f90
	$(FC) $(FFLAGS) -c gradient_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

gradient_gaussian_grid_saved.o: \
gradient_routines.o \
gradient_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c gradient_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_gradient_regular_grid.o: \
gradient_routines.o \
invert_gradient_regular_grid.f90
	$(FC) $(FFLAGS) -c invert_gradient_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_gradient_regular_grid_saved.o: \
gradient_routines.o \
invert_gradient_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_gradient_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_gradient_gaussian_grid.o: \
gradient_routines.o \
invert_gradient_gaussian_grid.f90
	$(FC) $(FFLAGS) -c invert_gradient_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_gradient_gaussian_grid_saved.o: \
gradient_routines.o \
invert_gradient_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_gradient_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

divergence_routines.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
vector_synthesis_routines.o \
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

invert_divergence_regular_grid.o: \
divergence_routines.o \
invert_divergence_regular_grid.f90
	$(FC) $(FFLAGS) -c invert_divergence_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_divergence_regular_grid_saved.o: \
divergence_routines.o \
invert_divergence_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_divergence_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_divergence_gaussian_grid.o: \
divergence_routines.o \
invert_divergence_gaussian_grid.f90
	$(FC) $(FFLAGS) -c invert_divergence_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_divergence_gaussian_grid_saved.o: \
divergence_routines.o \
invert_divergence_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_divergence_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vorticity_routines.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
vector_synthesis_routines.o \
vorticity_routines.f90
	$(FC) $(FFLAGS) -c vorticity_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vorticity_regular_grid.o: \
vorticity_routines.o \
vorticity_regular_grid.f90
	$(FC) $(FFLAGS) -c vorticity_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vorticity_regular_grid_saved.o: \
vorticity_routines.o \
vorticity_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c vorticity_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vorticity_gaussian_grid.o: \
vorticity_routines.o \
vorticity_gaussian_grid.f90
	$(FC) $(FFLAGS) -c vorticity_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vorticity_gaussian_grid_saved.o: \
vorticity_routines.o \
vorticity_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c vorticity_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vorticity_regular_grid.o: \
vorticity_routines.o \
invert_vorticity_regular_grid.f90
	$(FC) $(FFLAGS) -c invert_vorticity_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vorticity_regular_grid_saved.o: \
vorticity_routines.o \
invert_vorticity_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_vorticity_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vorticity_gaussian_grid.o: \
vorticity_routines.o \
invert_vorticity_gaussian_grid.f90
	$(FC) $(FFLAGS) -c invert_vorticity_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vorticity_gaussian_grid_saved.o: \
vorticity_routines.o \
invert_vorticity_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_vorticity_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_laplacian_routines.o: \
spherepack_precision.o \
scalar_synthesis_routines.o \
scalar_laplacian_routines.f90
	$(FC) $(FFLAGS) -c scalar_laplacian_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_laplacian_regular_grid.o: \
scalar_laplacian_routines.o \
scalar_laplacian_regular_grid.f90
	$(FC) $(FFLAGS) -c scalar_laplacian_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_laplacian_regular_grid_saved.o: \
scalar_laplacian_routines.o \
scalar_laplacian_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c scalar_laplacian_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_laplacian_gaussian_grid.o: \
scalar_laplacian_routines.o \
scalar_laplacian_gaussian_grid.f90
	$(FC) $(FFLAGS) -c scalar_laplacian_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

scalar_laplacian_gaussian_grid_saved.o: \
scalar_laplacian_routines.o \
scalar_laplacian_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c scalar_laplacian_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_scalar_laplacian_regular_grid.o: \
scalar_laplacian_routines.o \
invert_scalar_laplacian_regular_grid.f90
	$(FC) $(FFLAGS) -c invert_scalar_laplacian_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_scalar_laplacian_regular_grid_saved.o: \
scalar_laplacian_routines.o \
invert_scalar_laplacian_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_scalar_laplacian_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_scalar_laplacian_gaussian_grid.o: \
scalar_laplacian_routines.o \
invert_scalar_laplacian_gaussian_grid.f90
	$(FC) $(FFLAGS) -c invert_scalar_laplacian_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_scalar_laplacian_gaussian_grid_saved.o: \
scalar_laplacian_routines.o \
invert_scalar_laplacian_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_scalar_laplacian_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_laplacian_routines.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
vector_laplacian_routines.f90
	$(FC) $(FFLAGS) -c vector_laplacian_routines.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_laplacian_regular_grid.o: \
vector_laplacian_routines.o \
vector_laplacian_regular_grid.f90
	$(FC) $(FFLAGS) -c vector_laplacian_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_laplacian_regular_grid_saved.o: \
vector_laplacian_routines.o \
vector_laplacian_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c vector_laplacian_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_laplacian_gaussian_grid.o: \
vector_laplacian_routines.o \
vector_laplacian_gaussian_grid.f90
	$(FC) $(FFLAGS) -c vector_laplacian_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

vector_laplacian_gaussian_grid_saved.o: \
vector_laplacian_routines.o \
vector_laplacian_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c vector_laplacian_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vector_laplacian_regular_grid.o: \
vector_laplacian_routines.o \
invert_vector_laplacian_regular_grid.f90
	$(FC) $(FFLAGS) -c invert_vector_laplacian_regular_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vector_laplacian_regular_grid_saved.o: \
vector_laplacian_routines.o \
invert_vector_laplacian_regular_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_vector_laplacian_regular_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vector_laplacian_gaussian_grid.o: \
vector_laplacian_routines.o \
invert_vector_laplacian_gaussian_grid.f90
	$(FC) $(FFLAGS) -c invert_vector_laplacian_gaussian_grid.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

invert_vector_laplacian_gaussian_grid_saved.o: \
vector_laplacian_routines.o \
invert_vector_laplacian_gaussian_grid_saved.f90
	$(FC) $(FFLAGS) -c invert_vector_laplacian_gaussian_grid_saved.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtec.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
idvtec.f90
	$(FC) $(FFLAGS) -c idvtec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtes.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
idvtes.f90
	$(FC) $(FFLAGS) -c idvtes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtgc.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
idvtgc.f90
	$(FC) $(FFLAGS) -c idvtgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

idvtgs.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
idvtgs.f90
	$(FC) $(FFLAGS) -c idvtgs.f90 $(LINKS)\
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
vector_synthesis_routines.o \
isfvpec.f90
	$(FC) $(FFLAGS) -c isfvpec.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpes.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
isfvpes.f90
	$(FC) $(FFLAGS) -c isfvpes.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpgc.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
isfvpgc.f90
	$(FC) $(FFLAGS) -c isfvpgc.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

isfvpgs.o: \
spherepack_precision.o \
vector_synthesis_routines.o \
isfvpgs.f90
	$(FC) $(FFLAGS) -c isfvpgs.f90 $(LINKS)\
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

type_RealHarmonicCoefficients.o: \
spherepack_precision.o \
type_RealHarmonicCoefficients.f90
	$(FC) $(FFLAGS) -c type_RealHarmonicCoefficients.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_VectorHarmonicCoefficients.o: \
spherepack_precision.o \
type_RealHarmonicCoefficients.o \
type_VectorHarmonicCoefficients.f90
	$(FC) $(FFLAGS) -c type_VectorHarmonicCoefficients.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Workspace.o: \
spherepack_precision.o \
type_RealHarmonicCoefficients.o \
type_VectorHarmonicCoefficients.o \
type_Workspace.f90
	$(FC) $(FFLAGS) -c type_Workspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RegularWorkspace.o: \
spherepack_precision.o \
type_Workspace.o\
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vector_synthesis_routines.o \
type_RegularWorkspace.f90
	$(FC) $(FFLAGS) -c type_RegularWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_GaussianWorkspace.o: \
spherepack_precision.o \
type_Workspace.o\
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vector_synthesis_routines.o \
type_GaussianWorkspace.f90
	$(FC) $(FFLAGS) -c type_GaussianWorkspace.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Vector3D.o: \
spherepack_precision.o \
type_Vector3D.f90
	$(FC) $(FFLAGS) -c type_Vector3D.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_TrigonometricTable.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_TrigonometricTable.f90
	$(FC) $(FFLAGS) -c type_TrigonometricTable.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_SphericalUnitVectors.o: \
spherepack_precision.o \
type_SphericalGrid.o\
type_TrigonometricTable.o\
type_Vector3D.o \
type_SphericalUnitVectors.f90
	$(FC) $(FFLAGS) -c type_SphericalUnitVectors.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_Sphere.o: \
spherepack_precision.o \
type_TrigonometricTable.o \
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
type_TrigonometricTable.o \
type_Vector3D.o \
type_GaussianGrid.o\
type_GaussianWorkspace.o\
type_SphericalUnitVectors.o \
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vector_synthesis_routines.o \
type_GaussianSphere.f90
	$(FC) $(FFLAGS) -c type_GaussianSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

type_RegularSphere.o: \
spherepack_precision.o \
type_Sphere.o \
type_TrigonometricTable.o \
type_Vector3D.o \
type_RegularGrid.o\
type_RegularWorkspace.o\
type_SphericalUnitVectors.o \
scalar_analysis_routines.o \
scalar_synthesis_routines.o \
vector_analysis_routines.o \
vector_synthesis_routines.o \
type_RegularSphere.f90
	$(FC) $(FFLAGS) -c type_RegularSphere.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

spherepack.o: \
spherepack_precision.o \
divergence_routines.o \
gaussian_latitudes_and_weights_routines.o \
coordinate_transfer_routines.o \
gradient_routines.o \
idvtec.o \
idvtes.o \
idvtgc.o \
idvtgs.o \
icosahedral_geodesic_routines.o \
isfvpec.o \
isfvpes.o \
isfvpgc.o \
isfvpgs.o \
sfvpec.o \
sfvpes.o \
sfvpgc.o \
sfvpgs.o \
scalar_analysis_routines.o \
scalar_laplacian_routines.o \
scalar_projection_routines.o \
scalar_synthesis_routines.o \
sshifte.o \
trssph.o \
trvsph.o \
vector_analysis_routines.o \
vector_synthesis_routines.o \
vector_laplacian_routines.o \
vorticity_routines.o \
vshifte.o \
colatitudinal_derivative_routines.o \
type_RealPeriodicFastFourierTransform.o \
type_FFTpack.o \
type_LegendreAux.o \
type_Vector3D.o \
type_GaussianSphere.o \
type_RegularSphere.o \
spherepack.f90
	$(FC) $(FFLAGS) -c spherepack.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(MOD_DIRECTORY)

main.o: \
spherepack.o \
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
	rm -f $(OBJS) $(LIB) $(MOD_DIRECTORY)/*.mod $(MOD_DIRECTORY)/*.smod
	rm -f $(EXECUTABLE)