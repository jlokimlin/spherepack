include ../compiler_flags.inc
include ./relative_directories.inc

# Add module paths

spherepack_wrapperModulePath = ../../obj

modulePath = \
-I$(spherepack_wrapperModulePath) \
-I$(obj_directory)

# Add the required internal libraries
externalLibraryPath = /usr/local/lib

libInclude = \
-L$(externalLibraryPath) \
-lspherepack_wrapper -lspherepack \
$(modulePath)

executable = ./exe
	 
test.o: \
test.f90
	$(fcomp) -c test.f90 $(libInclude) -o $(obj_directory)/$@

exec: $(objects)
	$(fcomp) -o $(executable) $(objects_with_prefix) $(libInclude)

run:
	make clean
	make exec
	./exe
	 
clean:
	rm -f $(objects)
	rm -f $(executable)