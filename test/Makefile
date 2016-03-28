
include ../make.inc

SPHEREPACK_WRAPPER_MODULE_PATH = $(LIB_DIRECTORY)

LINKS = \
-L$(LIB_DIRECTORY) -lspherepack_wrapper \
-I$(SPHEREPACK_WRAPPER_MODULE_PATH) \
-I$(OBJ_DIRECTORY)

SRC=$(subst .o,.f90, $(OBJ))

EXES=$(subst .o,.exe, $(OBJ))

ENFORCE_64_BIT_PRECISION = \
-fdefault-real-8 \
-fdefault-double-8

F90 := gfortran $(ENFORCE_64_BIT_PRECISION)

run: $(EXES)

$(EXES): $(SRC)

%.exe : %.f90
	rm -f $@
	$(F90) $(FFLAGS) $< -o $@ $(LINKS)
	./$@

clean:
	rm -f $(OBJS) $(EXES)
