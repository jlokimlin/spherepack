
include ../make.inc

LIB_OLD_SPHEREPACK = $(LIB_DIRECTORY)/libold_spherepack.a

ENFORCE_64_BIT_PRECISION = \
-fdefault-real-8 \
-fdefault-double-8

F90 = gfortran 
FFLAG = -O3 $(ENFORCE_64_BIT_PRECISION)

$(LIB_OLD_SPHEREPACK) : $(OBJS)
	$(AR) -rv $@ $?

$(OBJ_DIRECTORY)/%.o : %.f90
	$(F90) $(FFLAG) -c $< -o $@

clean:
	rm -f $(OBJS) $(LIB_OLD_SPHEREPACK)
