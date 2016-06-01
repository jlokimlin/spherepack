!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK version 3.2                   *
!     *                                                               *
!     *       A Package of Fortran Subroutines and Programs           *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!
! ... file helmsph.f
!
!     this file contains a program for solving the Helmholtz
!     equation with constant 1.0 on a ten degree grid on the full sphere
!
! ... required spherepack files
!
!     islapec.f, shaec.f, shsec.f, sphcom.f, hrfft.f
!
! ... description
!
!     let theta be latitude and phi be east longitude in radians.
!     and let
!
!
!       x = cos(theta)*sin(phi)
!       y = cos(theta)*cos(phi)
!       z = sint(theta)
!
!     be the cartesian coordinates corresponding to theta and phi.
!     on the unit sphere.  The exact solution
!
!        ue(theta,phi) = (1.+x*y)*exp(z)
!
!     is used to set the right hand side and compute error.
!
!
program helmsph

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        Sphere, &
        Regularsphere, &
        GaussianSphere

    ! Explicit typing only
    implicit none

    !----------------------------------------------------------------------
    ! Dictionary
    !----------------------------------------------------------------------
    class (Sphere), pointer :: sphere_dat
    !----------------------------------------------------------------------

    ! Cast to gaussian case
    allocate( GaussianSphere :: sphere_dat )
    call test_helmholtz_inversion(sphere_dat)
    deallocate( sphere_dat )

    ! Cast to regular case
    allocate( RegularSphere :: sphere_dat )
    call test_helmholtz_inversion(sphere_dat)
    deallocate( sphere_dat )

contains


    subroutine test_helmholtz_inversion(sphere_type)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: sphere_type
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip), parameter        :: NLONS = 36
        integer (ip), parameter        :: NLATS = NLONS/2 + 1
        integer (ip)                   :: i, j !! Counters
        real (wp)                      :: exact_solution(NLATS, NLONS)
        real (wp)                      :: approximate_solution(NLATS, NLONS)
        real (wp)                      :: source_term(NLATS, NLONS)
        real (wp), parameter           :: HELMHOLTZ_CONSTANT = 1.0_wp
        character (len=:), allocatable :: error_previous_platform
        !----------------------------------------------------------------------

        !
        !==> Set up workspace arrays
        !
        select type (sphere_type)
            type is (GaussianSphere)

            !  Initialize gaussian sphere object
            sphere_type = GaussianSphere(NLATS,NLONS)

            ! Allocate known error from previous platform
            allocate( error_previous_platform, source='     discretization error = 2.325553e-14' )

            type is (RegularSphere)

            ! Initialize regular sphere
            sphere_type = RegularSphere(NLATS, NLONS)

            ! Allocate known error from previous platform
            allocate( error_previous_platform, source='     discretization error = 1.202313e-14' )
        end select

        !
        !==> Set right hand side as helmholtz operator
        !    applied to ue = (1.+x*y)*exp(z)
        !
        associate( &
            ue => exact_solution, &
            rhs => source_term, &
            radial => sphere_type%unit_vectors%radial &
            )
            do j=1,NLONS
                do i=1,NLATS
                    associate( &
                        x => radial(i,j)%x, &
                        y => radial(i,j)%y, &
                        z => radial(i,j)%z &
                        )
                        ue(i,j) = (1.0_wp + x * y) * exp(z)
                        rhs(i,j) = -(x*y*(z*z+6.0_wp*(z+1.0_wp))+z*(z+2.0_wp))*exp(z)
                    end associate
                end do
            end do
        end associate


        associate( &
            xlmbda => HELMHOLTZ_CONSTANT, &
            rhs => source_term, &
            u => approximate_solution &
            )
            !
            !==> Solve Helmholtz equation on the sphere
            !
            call sphere_type%invert_helmholtz( xlmbda, rhs, u)

        end associate

        !
        !==> Compare ue with u
        !
        associate( &
            u => approximate_solution, &
            ue => exact_solution &
            )
            associate( err2 => norm2(u-ue) )
                !
                !==> Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     helmsph *** TEST RUN *** '
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(A)') '     Helmholtz approximation on a ten degree grid'
                write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
                write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(A)') error_previous_platform
                write( stdout, '(A)') '     The output from your computer is: '
                write( stdout, '(A,1pe15.6)') '     discretization error = ', err2
                write( stdout, '(A)' ) ''
            end associate
        end associate
        !
        !==> Release memory
        !
        call sphere_type%destroy()
        deallocate( error_previous_platform )

    end subroutine test_helmholtz_inversion


end program helmsph
