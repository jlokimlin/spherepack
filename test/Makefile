
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

FC = gfortran

FCFLAGS = -O3 $(ENFORCE_64_BIT_PRECISION)

run: $(EXES)

$(EXES): $(SRC)

%.exe : %.f90
	rm -f $@
	$(FC) $(FCFLAGS) $< -o $@ $(LINKS)
	./$@

clean:
	rm -f $(OBJS) $(EXES)
