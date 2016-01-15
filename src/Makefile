include ./compiler_flags.inc
include ./directories.inc

MODULE_PATH = $(OBJ_DIRECTORY)

EXTERNAL_LIBRARY_PATH = /usr/local/lib

LINKS = \
-lspherepack -L$(EXTERNAL_LIBRARY_PATH) \
-I$(MODULE_PATH)

EXECUTABLE = ./exe

type_vector_mod.o: \
type_vector_mod.f90
	$(FC) $(FFLAGS) -c type_vector_mod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

type_grid_mod.o: \
type_grid_mod.f90
	$(FC) $(FFLAGS) -c type_grid_mod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

type_workspace_mod.o: \
type_workspace_mod.f90
	$(FC) $(FFLAGS) -c type_workspace_mod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

type_sphere_mod.o: \
type_sphere_mod.f90
	$(FC) $(FFLAGS) -c type_sphere_mod.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@ -J$(OBJ_DIRECTORY)

spherepack_wrapper_mod.o: \
type_vector_mod.o \
type_grid_mod.o\
type_workspace_mod.o\
type_sphere_mod.o \
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
	ar rc libspherepack_wrapper.a \
	$(OBJ_DIRECTORY)/type_vector_mod.o \
	$(OBJ_DIRECTORY)/type_grid_mod.o \
	$(OBJ_DIRECTORY)/type_workspace_mod.o \
	$(OBJ_DIRECTORY)/type_sphere_mod.o \
	$(OBJ_DIRECTORY)/spherepack_wrapper_mod.o
	
	cp libspherepack_wrapper.a /usr/local/lib
	mv libspherepack_wrapper.a $(LIB_DIRECTORY)

all: 
	make run
	cd ./test; make run

run:
	make clean
	make exec
	make lib
	./exe
	 
clean:
	rm -f $(OBJ_DIRECTORY)/*.o $(OBJ_DIRECTORY)/*.mod
	rm -f $(EXECUTABLE)
