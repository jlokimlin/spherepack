
include ../make.inc

LIB_SPHEREPACK=../lib/libold_spherepack.a

ENFORCE_64_BIT_PRECISION = \
-fdefault-real-8 \
-fdefault-double-8

F90 := gfortran -O3 $(ENFORCE_64_BIT_PRECISION)

$(LIB_SPHEREPACK) : $(OBJS)
	$(AR) -rv $@ $?

$(OBJ_DIRECTORY)/%.o : %.f90
	$(F90) -c $< -o $@

clean:
	rm -f $(OBJS) 
