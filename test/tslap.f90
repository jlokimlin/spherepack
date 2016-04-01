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

    call test_scalar_laplacian_routines( gaussian_sphere )
    call test_scalar_laplacian_routines( regular_sphere )


contains

    subroutine test_scalar_laplacian_routines( sphere_type )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: sphere_type
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip), parameter        :: NLATS = 15
        integer (ip), parameter        :: NLONS = 22
        integer (ip), parameter        :: NSYNTHS = 3
        integer (ip)                   :: i, j, k !! Counters
        real (wp)                      :: scalar_function(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: approximate_laplacian(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: exact_laplacian(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: err2, se
        character (len=:), allocatable :: laplacian_error, inversion_error
        !----------------------------------------------------------------------

        !
        !==> Set up workspace arrays
        !
        select type (sphere_type)
            !
            !==> For gaussian sphere
            !
            class is (GaussianSphere)
            !  Initialize gaussian sphere object
            call sphere_type%create(nlat=NLATS, nlon=NLONS)
            ! Allocate known error from previous platform
            allocate( laplacian_error, source='     discretization error = 1.682230e-13' )
            allocate( inversion_error, source='     discretization error = 9.988022e-16' )
            !
            !==> For regular sphere
            !
            class is (RegularSphere)
            ! Initialize regular sphere
            call sphere_type%create(nlat=NLATS, nlon=NLONS)
            ! Allocate known error from previous platform
            allocate( laplacian_error, source='     discretization error = 1.450713e-13' )
            allocate( inversion_error, source='     discretization error = 9.666880e-16' )
        end select

        !
        !==> test all analysis and synthesis subroutines
        !    set scalar field as polynomial in x,y,z
        !    restricted to the sphere
        !
        associate( &
            s => scalar_function, &
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
                                    s(i,j,k) = x + y
                                    sclpe(i,j,k) = -2.0_wp * (x + y)
                                case(2)
                                    s(i,j,k) = x+z
                                    sclpe(i,j,k) = -2.0_wp * (x + z)
                                case(3)
                                    s(i,j,k) = y+z
                                    sclpe(i,j,k) = -2.0_wp * (y + z)
                            end select
                        end associate
                    end do
                end do
            end do
        end associate

        !
        !==> Compute scalar laplacian
        !
        do k = 1, NSYNTHS
            associate( &
                s => scalar_function(:,:,k), &
                sclp => approximate_laplacian(:,:,k) &
                )
                call sphere_type%get_laplacian(s, sclp)
            end associate
        end do

        !
        !==> Compute discretization error in sclp
        !
        err2 = 0.0_wp
        associate( &
            sclp => approximate_laplacian, &
            sclpe => exact_laplacian &
            )
            do k = 1, NSYNTHS
                do j = 1, NLONS
                    do i = 1, NLATS
                        err2 = err2 + (sclpe(i,j,k)-sclp(i,j,k))**2
                    end do
                end do
            end do
            err2 = sqrt(err2/size(sclp))
        end associate

        !
        !==> Print earlier output from platform with 64-bit floating point
        !    arithmetic followed by the output from this computer
        !
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     tslap *** TEST RUN *** '
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
        write( stdout, '(A)') '     scalar laplacian approximation'
        write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') laplacian_error
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,1pe15.6)') '     discretization error = ', err2
        write( stdout, '(A)' ) ''

        !
        !==> invert scalar laplacian
        !
        do k = 1, NSYNTHS
            associate( &
                sclpe => exact_laplacian(:,:,k), &
                s => scalar_function(:,:,k) &
                )
                call sphere_type%invert_laplacian(sclpe, s)
            end associate
        end do

        !
        !==> Compare s with original
        !
        err2 = 0.0
        associate( &
            s => scalar_function, &
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
                                case (1)
                                    se = x+y
                                case (2)
                                    se = x+z
                                case (3)
                                    se = y+z
                            end select
                            err2 = err2+(s(i,j,k) - se)**2
                        end associate
                    end do
                end do
            end do
            err2 = sqrt(err2/size(s))
        end associate

        !
        !==> Print earlier output from platform with 64-bit floating point
        !    arithmetic followed by the output from this computer
        !
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     tslap *** TEST RUN *** '
        write( stdout, '(A)') ''
        write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
        write( stdout, '(A)') '     scalar laplacian inversion'
        write( stdout, '(2(A,I2))') '     nlat = ', NLATS,' nlon = ', NLONS
        write( stdout, '(A)') '     Previous 64 bit floating point arithmetic result '
        write( stdout, '(A)') laplacian_error
        write( stdout, '(A)') '     The output from your computer is: '
        write( stdout, '(A,1pe15.6)') '     discretization error = ', err2
        write( stdout, '(A)' ) ''

        !
        !==> Release memory
        !
        call sphere_type%destroy()
        deallocate( laplacian_error )


    end subroutine test_scalar_laplacian_routines

end program tslap


