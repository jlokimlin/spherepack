# **spherepack\_wrapper**

This Fortran project is a modernization of NCAR's SPHEREPACK3.2 library. The original work, written in FORTRAN 77, was heavily refactored to incorporate features of modern Fortran (2008+). The arduous initialization procedures for analysis and synthesis are now confined to the polymorphic class variable **Sphere**. This object-oriented approach hides the various workspace arrays from the user.

This project is still a work in progress and mainly works with gaussian grids; regular (equally-spaced) grids are only partially supported.
 

-----------------------------------------------------------------------------

## Requirements
* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------


## To build the project

Type the following command line arguments
```
git clone https://github.com/jlokimlin/spherepack_wrapper.git

cd spherepack_wrapper; make all
```

-----------------------------------------------------------------------------

## TODO
* Replace old-style **do** loops
* Replace all instances of **go to** statements with **exit**, **cycle** and **select case**
* Introduce parameterized kinds to remove the compiler flags **-fdefault-real-8 -fdefault-double-8**