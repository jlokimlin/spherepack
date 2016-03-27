
include ../make.inc

LIB_SPHEREPACK=../lib/libspherepack.a

ENFORCE_64_BIT_PRECISION = \
-fdefault-real-8 \
-fdefault-double-8

F90 := gfortran $(ENFORCE_64_BIT_PRECISION)

$(LIB_SPHEREPACK) : $(OBJS)
	$(AR) -rv $@ $?

$(OBJ_DIRECTORY)/%.o : %.f90
	$(F90) -c $< -o $@

clean:
	rm -f $(OBJS) 
