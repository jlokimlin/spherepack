program helmsph

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use modern_spherepack_library, only: &
        Sphere, &
        Regularsphere, &
        GaussianSphere

    ! Explicit typing only
    implicit none

    !----------------------------------------------------------------------
    ! Dictionary
    !----------------------------------------------------------------------
    type (GaussianSphere) :: gaussian_sphere
    type (RegularSphere)  :: regular_sphere
    !----------------------------------------------------------------------

    call test_helmsph( gaussian_sphere )
    call test_helmsph( regular_sphere )


contains


    subroutine test_helmsph( sphere_type )
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
        !     *       A Package of Fortran77 Subroutines and Programs         *
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
        real (wp)                      :: approximate_solution(NLATS, NLONS)
        real (wp)                      :: source_term(NLATS, NLONS)
        real (wp)                      :: helmholtz_constant, discretization_error
        character (len=:), allocatable :: prev_error
        !----------------------------------------------------------------------

        !
        !==> Set up workspace arrays
        !
        select type (sphere_type)
            class is (GaussianSphere)
            ! Create gaussian sphere
            call sphere_type%create(nlat=NLATS, nlon=NLONS, isym=0, isynt=1)
            ! Allocate prev known
            allocate( prev_error, source='     discretization error = 3.552714e-15' )
            class is (RegularSphere)
            ! Create regular sphere
            call sphere_type%create(nlat=NLATS, nlon=NLONS, isym=0, isynt=1)
            ! Allocate prev known
            allocate( prev_error, source='     discretization error = 2.331468e-15' )
        end select

        !
        !==> Set helmholtz constant
        !
        helmholtz_constant = 1.0_wp

        ! Set right hand side as helmholtz operator
        ! applied to ue = (1.+x*y)*exp(z)
        associate( &
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
                        rhs(i,j) = -(x*y*(z*z+6.0_wp*(z+1.0_wp))+z*(z+2.0_wp))*exp(z)
                    end associate
                end do
            end do
        end associate

        !
        !==> Solve Helmholtz equation on the sphere in u
        !
        associate( &
            xlmbda => helmholtz_constant, &
            rhs => source_term, &
            u => approximate_solution &
            )
            call sphere_type%invert_helmholtz( xlmbda, rhs, u)
        end associate

        !
        !==> Compute and print maximum error
        !
        associate( &
            err_max => discretization_error, &
            u => approximate_solution, &
            radial => sphere_type%unit_vectors%radial &
            )
            ! Initialize error
            err_max = 0.0_wp
            do j=1,NLONS
                do i=1,NLATS
                    associate( &
                        ! Associate radial components
                        x => radial(i,j)%x, &
                        y => radial(i,j)%y, &
                        z => radial(i,j)%z &
                        )
                        ! Set exact solution
                        associate( ue => (1.0_wp + x * y) * exp(z) )
                            err_max = max(err_max,abs(u(i,j)-ue))
                        end associate
                    end associate
                end do
            end do
        end associate

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
        write( stdout, '(A)') prev_error
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,1pe15.6)') '     discretization error = ', &
            discretization_error
        write( stdout, '(A)' ) ''

        !
        !==> Release memory
        !
        call sphere_type%destroy()
        deallocate( prev_error )

    end subroutine test_helmsph


end program helmsph

