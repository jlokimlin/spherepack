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
!
!     a program for testing all scalar analysis and synthesis subroutines
!
program tsha

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
    type (GaussianSphere) :: gaussian_sphere
    type (RegularSphere)  :: regular_sphere
    !----------------------------------------------------------------------

    call test_analysis_and_synthesis_routines( gaussian_sphere )
    call test_analysis_and_synthesis_routines( regular_sphere )


contains

    subroutine test_analysis_and_synthesis_routines( sphere_type )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: sphere_type
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip), parameter        :: NLONS = 128
        integer (ip), parameter        :: NLATS = NLONS/2 + 1
        integer (ip), parameter        :: NSYNTHS = 3
        integer (ip)                   :: i, j, k !! Counters
        real (wp)                      :: original_scalar_function(NLATS,NLONS,NSYNTHS)
        real (wp)                      :: approximate_scalar_function(NLATS,NLONS,NSYNTHS)
        character (len=:), allocatable :: error_previous_platform
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
            allocate( error_previous_platform, source='     discretization error = 3.375078e-14' )
            !
            !==> For regular sphere
            !
            class is (RegularSphere)

            ! Initialize regular sphere
            call sphere_type%create(nlat=NLATS, nlon=NLONS)

            ! Allocate known error from previous platform
            allocate( error_previous_platform, source='     discretization error = 2.664535e-14' )
        end select

        !
        !==> Test all analysis and synthesis subroutines.
        !    Set scalar field as (x*y*z)**k) restricted to the sphere
        !
        associate( &
            se => original_scalar_function, &
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
                                    se(i,j,k) = exp( x + y + z )
                                case default
                                    se(i,j,k) = (x*y*z)**k
                            end select
                        end associate
                    end do
                end do
            end do
        end associate

        !
        !==> Perform analysis then synthesis
        !
        do k = 1, NSYNTHS
            associate( &
                se => original_scalar_function(:,:,k), &
                s => approximate_scalar_function(:,:,k) &
                )

                ! Analyse function into (real) coefficients
                call sphere_type%perform_complex_analysis(se)

                ! Synthesize function from (real) coefficients
                call sphere_type%perform_complex_synthesis(s)
            end associate
        end do
        !
        !==> Compute discretization error
        !
        associate( &
            se => original_scalar_function, &
            s => approximate_scalar_function &
            )
            associate( err2 => maxval(abs(s - se)) )
                !
                !==> Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                !
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     tsha *** TEST RUN *** '
                write( stdout, '(A)') ''
                write( stdout, '(A)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(A)') '     Testing scalar analysis and synthesis'
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

    end subroutine test_analysis_and_synthesis_routines

end program tsha
