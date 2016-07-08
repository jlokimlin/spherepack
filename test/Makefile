
include ../make.inc

LINKS = \
-L$(LIB_DIRECTORY) -l$(LIB_NAME) \
-I$(MOD_DIRECTORY) \
-I$(OBJ_DIRECTORY)

SRC=$(subst .o,.f90, $(OBJ))

EXES=$(subst .o,.exe, $(OBJ))

ENFORCE_64_BIT_PRECISION = \
-fdefault-real-8 \
-fdefault-double-8

F90 = gfortran

FFLAGS = -O3 $(ENFORCE_64_BIT_PRECISION)

all : clean run

run: $(EXES)

$(EXES): $(SRC)

%.exe : %.f90
	rm -f $@
	$(F90) $(FFLAGS) $< -o $@ $(LINKS)
	./$@

clean:
	rm -f $(OBJS) $(EXES) ./*.mod
