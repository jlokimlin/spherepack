include ../fcomp.inc
include ./relative_directories.inc
include ../libInclude.inc

# Add the required internal libraries

libInclude += \
-I$(obj_directory) \
-lspherepack_wrapper -L$(externalLibraryPath)

executable = ./exe
	 
test.o: \
test.f90
	$(fcomp) -c test.f90 $(libInclude)\
	 -o $(obj_directory)/$@

exec: $(objects)
	$(fcomp) -o $(executable) $(objects_with_prefix) $(libInclude)

run:
	make clean
	make exec
	./exe
	 
clean:
	rm -f $(objects)
	rm -f $(executable)