include ./fcomp.inc
include ./directories.inc

# Add the required internal libraries
externalLibraryPath = /usr/local/lib

libInclude = \
-lspherepack -L$(externalLibraryPath) \
-I$(obj_directory)

executable = ./exe

type_vector_mod.o: \
type_vector_mod.f90
	$(fcomp) -c type_vector_mod.f90 $(libInclude)\
	 -o $(obj_directory)/$@ -J$(obj_directory)

type_grid_mod.o: \
type_grid_mod.f90
	$(fcomp) -c type_grid_mod.f90 $(libInclude)\
	 -o $(obj_directory)/$@ -J$(obj_directory)

type_workspace_mod.o: \
type_workspace_mod.f90
	$(fcomp) -c type_workspace_mod.f90 $(libInclude)\
	 -o $(obj_directory)/$@ -J$(obj_directory)

type_sphere_mod.o: \
type_sphere_mod.f90
	$(fcomp) -c type_sphere_mod.f90 $(libInclude)\
	 -o $(obj_directory)/$@ -J$(obj_directory)

spherepack_wrapper_mod.o: \
type_vector_mod.o \
type_grid_mod.o\
type_workspace_mod.o\
type_sphere_mod.o \
spherepack_wrapper_mod.f90
	$(fcomp) -c spherepack_wrapper_mod.f90 $(libInclude)\
	 -o $(obj_directory)/$@ -J$(obj_directory)

main.o: \
spherepack_wrapper_mod.o \
main.f90
	$(fcomp) -c main.f90 $(libInclude)\
	 -o $(obj_directory)/$@

exec: $(objects)
	$(fcomp) -o $(executable) $(objects_with_prefix) $(libInclude)

build_lib: 
	ar rc libspherepack_wrapper.a \
	$(obj_directory)/type_vector_mod.o \
	$(obj_directory)/type_grid_mod.o \
	$(obj_directory)/type_workspace_mod.o \
	$(obj_directory)/type_sphere_mod.o \
	$(obj_directory)/spherepack_wrapper_mod.o
	
	cp libspherepack_wrapper.a /usr/local/lib
	mv libspherepack_wrapper.a $(lib_directory)

run:
	make clean
	make exec
	make build_lib
	./exe
	 
clean:
	rm -f $(objects)
	rm -f $(executable)
