# **spherepack\_wrapper**

This Fortran project implements a 64-bit precision object-oriented wrapper for NCAR's SPHEREPACK3.2.


https://www2.cisl.ucar.edu/resources/legacy/spherepack


This project is still a work in progress and mainly works with gaussian grids; regular (equally-spaced) grids are only partially supported.


-----------------------------------------------------------------------------

## Requirements
* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran
* NCAR's SPHEREPACK3.2 https://github.com/jlokimlin/spherepack3.2.git

-----------------------------------------------------------------------------


## To build the project

Type the following command line arguments
```
git clone https://github.com/jlokimlin/spherepack_wrapper.git

cd spherepack_wrapper; make all
```