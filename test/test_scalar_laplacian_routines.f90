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
!     a program for testing the scalar laplacian routines
!
!     (1) set a scalar field s as poly in x,y,z restricted to sphere
!
!     (2) compute scalar laplacian in array sclp
!
!     (3) compare (2) with analytic scalar laplacian in sclpe
!
!     (4) compute the inverse  of (2) and compare with (1)
!
program tslap

    use, intrinsic :: ISO_Fortran_env, only: &
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
    class(Sphere), allocatable :: sphere_dat
    !----------------------------------------------------------------------

    !
    !  Test gaussian case
    !
    allocate( GaussianSphere :: sphere_dat )

    call test_scalar_laplacian_routines(sphere_dat)

    deallocate( sphere_dat )

    !
    !  Test regular case
    !
    allocate( RegularSphere :: sphere_dat )

    call test_scalar_laplacian_routines(sphere_dat)

    deallocate( sphere_dat )



contains



    subroutine test_scalar_laplacian_routines(sphere_type)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: sphere_type
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip), parameter        :: NLATS = 15
        integer(ip), parameter        :: NLONS = 22
        integer(ip), parameter        :: NSYNTHS = 3
        integer(ip)                   :: i, j, k !! Counters
        real(wp)                      :: original_scalar_function(NLATS,NLONS,NSYNTHS)
        real(wp)                      :: exact_laplacian(NLATS,NLONS,NSYNTHS)
        real(wp)                      :: approximate_scalar_function(NLATS,NLONS,NSYNTHS)
        real(wp)                      :: approximate_laplacian(NLATS,NLONS,NSYNTHS)
        character(len=:), allocatable :: laplacian_error, inversion_error
        !----------------------------------------------------------------------

        !
        !  Set up workspace arrays
        !
        select type(sphere_type)
            type is (GaussianSphere)

            !  Initialize gaussian sphere object
            sphere_type = GaussianSphere(NLATS, NLONS)

            ! Allocate known error from previous platform
            allocate( laplacian_error, source='     discretization error = 1.953993e-13' )
            allocate( inversion_error, source='     discretization error = 6.661338e-16' )

            type is (RegularSphere)

            ! Initialize regular sphere
            sphere_type = RegularSphere(NLATS, NLONS)

            ! Allocate known error from previous platform
            allocate( laplacian_error, source='     discretization error = 1.003642e-13' )
            allocate( inversion_error, source='     discretization error = 8.881784e-16' )
        end select

        !
        !  test all analysis and synthesis subroutines
        !    set scalar field as polynomial in x,y,z
        !    restricted to the sphere
        !
        associate( &
            se => original_scalar_function, &
            sclpe => exact_laplacian, &
            radial => sphere_type%unit_vectors%radial &
            )
            do k=1,NSYNTHS
                do j=1,NLONS
                    do i=1,NLATS
                        associate( &
                            x => radial(i,j)%x, &
                            y => radial(i,j)%y, &
                            z => radial(i,j)%z &
                            )
                            select case (k)
                                case(1)
                                    se(i,j,k) = x + y
                                    sclpe(i,j,k) = -2.0_wp * (x + y)
                                case(2)
                                    se(i,j,k) = x + z
                                    sclpe(i,j,k) = -2.0_wp * (x + z)
                                case(3)
                                    se(i,j,k) = y + z
                                    sclpe(i,j,k) = -2.0_wp * (y + z)
                            end select
                        end associate
                    end do
                end do
            end do
        end associate

        !
        !  Compute scalar laplacian
        !
        do k = 1, NSYNTHS
            associate( &
                se => original_scalar_function(:,:,k), &
                sclp => approximate_laplacian(:,:,k) &
                )
                call sphere_type%get_laplacian(se, sclp)
            end associate
        end do

        !
        !  Compute discretization error in sclp
        !
        associate( &
            sclp => approximate_laplacian, &
            sclpe => exact_laplacian &
            )
            associate( err2 => maxval(abs(sclp - sclpe)) )
                !
                !  Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     tslap *** TEST RUN *** '
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(a)') '     scalar laplacian approximation'
                write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
                write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(a)') laplacian_error
                write( stdout, '(a)') '     The output from your computer is: '
                write( stdout, '(a,1pe15.6)') '     discretization error = ', err2
                write( stdout, '(a)' ) ''
            end associate
        end associate
        !
        !  invert scalar laplacian
        !
        do k = 1, NSYNTHS
            associate( &
                sclpe => exact_laplacian(:,:,k), &
                s => approximate_scalar_function(:,:,k) &
                )
                call sphere_type%invert_laplacian(sclpe, s)
            end associate
        end do

        !
        !  Compare s with se
        !
        associate( &
            s => approximate_scalar_function, &
            se => original_scalar_function &
            )
            associate( err2 => maxval(abs(s - se)) )
                !
                !  Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     tslap *** TEST RUN *** '
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(a)') '     scalar laplacian inversion'
                write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
                write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(a)') inversion_error
                write( stdout, '(a)') '     The output from your computer is: '
                write( stdout, '(a,1pe15.6)') '     discretization error = ', err2
                write( stdout, '(a)' ) ''
            end associate
        end associate
        !
        !  Release memory
        !
        deallocate( laplacian_error )
        deallocate( inversion_error )
        call sphere_type%destroy()

    end subroutine test_scalar_laplacian_routines

end program tslap


