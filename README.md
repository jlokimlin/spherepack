fortran95 module interface to spherepack 
(http://www.scd.ucar.edu/css/software/spherepack)

works with both gaussian and regularly spaced grids.

See module source for documentation.

Prerequisites:  Spherepack library must be installed.  Build
single and double precision versions of lib using '-r8" fortran compiler
flag (or whatever option on your compiler automatically promotes real*4
to real*8).

To compile and test:

1) edit Makefile and change compiler name and flags, and path to find
spherepack libs.
2) make; make test compiles module and runs test program in single precision.
3) make dp; make testdp runs in double precision (only if your compiler
has a flag to automagically promote real to real*8).


-Jeff Whitaker
Jeffrey.S.Whitaker@noaa.gov
Sept, 2002