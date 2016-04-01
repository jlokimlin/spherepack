# **modern\_spherepack**

This Fortran project is a object-oriented programming (OOP) modernization of NCAR's SPHEREPACK3.2 library. The original work, written in FORTRAN 77, was heavily refactored to incorporate features of modern Fortran (2008+). 


The arduous initialization procedures for analysis and synthesis are now confined to the polymorphic class variable **Sphere**. The OOP approach hides the various workspace arrays from the user.

This project is still a work in progress and mainly works with gaussian grids; regular (equally-spaced) grids are only partially supported.
 

-----------------------------------------------------------------------------

## What is spherepack?

A collection of Fortran programs for computing common spherical differential operators including divergence, vorticity, latitudinal derivatives, gradients, the Laplacian of both scalar and vector functions, and the inverses of these operators.

For example, given divergence and vorticity, the package can be used to compute velocity components, then the Laplacian inverse can be used to solve the scalar and vector Poisson equations. The package also contains routines for computing the associated Legendre functions, Gaussian points and weights, multiple fast Fourier transforms, and for converting scalar and vector fields between geophysical and mathematical spherical coordinates.

Test programs are provided for solving these equations. Each program serves two purposes: as a template to guide you in writing your own codes utilizing the spherepack routines, and as a demonstration on your computer that you can correctly produce spherepack executables.

-----------------------------------------------------------------------------

## Usage

```fortran

    use spherepack_wrapper_library, only: &
        GaussianSphere

    ! Explicit typing only
    implicit none
    
    type (GaussianSphere)  :: sphere_dat
    real (wp), allocatable :: scalar_function(:,:)
    real (wp), allocatable :: laplacian(:,:) 
    real (wp), allocatable :: solution(:,:)
    
    ! Create object
    call sphere_dat%create(nlat=19, nlon=36)
    
    !.... generate some data
    
    ! Compute laplacian on sphere
    call sphere_dat%get_laplacian( scalar_function, laplacian )
    
    ! Invert laplacian on sphere
    call sphere_dat%invert_laplacian( laplacian, solution )
    
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
```
git clone https://github.com/jlokimlin/modern_spherepack.git

cd modern_spherepack; make all
```

-----------------------------------------------------------------------------


## TODO
* Replace old-style **do** loops
* Replace all instances of **go to** statements with **exit**, **cycle** and **select case**
* Introduce parameterized kinds to remove the compiler flags **-fdefault-real-8 -fdefault-double-8**