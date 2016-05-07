# **spherepack - A modern Fortran (2008+) library of spherical harmonic transforms**

An object-oriented programming (OOP) modernization of NCAR's SPHEREPACK3.2. 

The original work, written in FORTRAN 77, was heavily refactored to incorporate features of modern Fortran (2008+). The arduous initialization procedures for analysis and synthesis are now confined to the polymorphic class variable **Sphere**. The OOP approach hides the various workspace arrays from the user.

-----------------------------------------------------------------------------

## What is spherepack?

A collection of Fortran programs for computing common spherical differential operators including divergence, vorticity, latitudinal derivatives, gradients, the Laplacian of both scalar and vector functions, and the inverses of these operators.

For example, given divergence and vorticity, the package can be used to compute velocity components, then the Laplacian inverse can be used to solve the scalar and vector Poisson equations. The package also contains routines for computing the associated Legendre functions, Gaussian points and weights, multiple fast Fourier transforms, and for converting scalar and vector fields between geophysical and mathematical spherical coordinates.

Test programs are provided for solving these equations. Each program serves two purposes: as a template to guide you in writing your own codes utilizing the spherepack routines, and as a demonstration on your computer that you can correctly produce spherepack executables.

-----------------------------------------------------------------------------

## Usage

```fortran

    use spherepack_library, only: &
        GaussianSphere

    ! Explicit typing only
    implicit none
    
    type (GaussianSphere)  :: sphere_dat
    real (wp), allocatable :: scalar_function(:,:)
    real (wp), allocatable :: laplacian(:,:)
    real (wp), allocatable :: solution(:,:)
    
    ! Initialize object
    call sphere_dat%create(nlat=19, nlon=36)
    
    !.... generate some data
    
    ! Compute laplacian on sphere
    call sphere_dat%get_laplacian(scalar_function, laplacian)
    
    ! Invert laplacian on sphere
    call sphere_dat%invert_laplacian(laplacian, solution)
    
    ! Release memory
    call sphere_dat%destroy()

```

-----------------------------------------------------------------------------

## Requirements

* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------


## To build the project

Type the following command line arguments

```bash

	git clone https://github.com/jlokimlin/spherepack.git
	
	cd spherepack; make all
```

-----------------------------------------------------------------------------

## Contributing

This project is still a work in progress and anyone is free to contribute. 

-----------------------------------------------------------------------------

## TODO
* Enclose **assignment**(=) and **operator**(\*) inside **type**(ThreeDimensionalVector)
* Replace old-style **do** loops
* Replace all occurences of **go to** with **if-then-else**, **exit**, **cycle** and **select case**
* Universal instances of parameterized kinds *INT32* and *REAL64* to remove the compiler flags **-fdefault-real-8 -fdefault-double-8**
* Encapsulate all functions and subroutines inside modules.
* Improve documentation and build tools. 

-----------------------------------------------------------------------------


## Bibliography

[1] Swarztrauber, Paul N. "On computing the points and weights for Gauss--Legendre quadrature." *SIAM Journal on Scientific Computing* 24.3 (2003): 945-954.

[2] Swarztrauber, Paul N., and William F. Spotz. "Generalized discrete spherical harmonic transforms." *Journal of Computational Physics* 159.2 (2000): 213-230.

[3] Adams, John C., and Paul N. Swarztrauber. "SPHEREPACK 3.0: A model development facility." *Monthly Weather Review* 127.8 (1999): 1872-1878.

[4] Swarztrauber, Paul N. "Spectral transform methods for solving the shallow-water equations on the sphere." *Monthly Weather Review* 124.4 (1996): 730-744.

[5] Williamson, David L., et al. "A standard test set for numerical approximations to the shallow water equations in spherical geometry." *Journal of Computational Physics* 102.1 (1992): 211-224.



