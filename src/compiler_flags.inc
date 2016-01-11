# Set gfortran compiler options

fcomp = \
gfortran \
-O3 \
$(languageOptions) \
$(errorWarningOptions) \
$(codeGenerationOptions)

languageOptions = \
-fimplicit-none \
-fdefault-real-8 \
-fdefault-double-8

errorWarningOptions = \
-Wall

codeGenerationOptions = \
-fcheck=all \
-fbounds-check \
-fno-backtrace \
-fno-unsafe-math-optimizations \
-frounding-math \
-fsignaling-nans

# References:
# https://gcc.gnu.org/onlinedocs/gfortran/Code-Gen-Options.html#Code-Gen-Options
# https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html#Error-and-Warning-Options
# https://gcc.gnu.org/onlinedocs/gfortran/Option-Summary.html
