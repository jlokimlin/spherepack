!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK                               *
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
!
!     a program for testing all scalar analysis and synthesis subroutines
!
program test_analysis_and_synthesis_routines

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack, only: &
        wp, & ! Working precison
        ip, & ! Integer precision
        Sphere, &
        Regularsphere, &
        GaussianSphere

    ! Explicit typing only
    implicit none

    ! Dictionary
    class(Sphere), allocatable :: solver

    ! Test gaussian grid
    allocate (GaussianSphere :: solver)
    call test_case(solver)
    deallocate (solver)

    ! Test regular grid
    allocate (RegularSphere :: solver)
    call test_case(solver)
    deallocate (solver)

contains

    subroutine test_case(sphere_type)

        ! Dummy arguments
        class(Sphere), intent(inout)  :: sphere_type

        ! Local variables
        integer(ip), parameter        :: NLONS = 128
        integer(ip), parameter        :: NLATS = NLONS/2 + 1
        integer(ip), parameter        :: NSYNTHS = 3
        integer(ip)                   :: i, j, k ! Counters
        real(wp)                      :: original_scalar_function(NLATS, NLONS, NSYNTHS)
        real(wp)                      :: approximate_scalar_function(NLATS, NLONS, NSYNTHS)
        character(len=:), allocatable :: error_previous_platform

        !  Set up workspace arrays
        select type(sphere_type)
            type is (GaussianSphere)

            !  Initialize gaussian sphere object
            sphere_type = GaussianSphere(NLATS, NLONS)

            ! Allocate known error from previous platform
            allocate (error_previous_platform, &
                source='     discretization error = 3.375078e-14')

            type is (RegularSphere)

            ! Initialize regular sphere
            sphere_type = RegularSphere(NLATS, NLONS)

            ! Allocate known error from previous platform
            allocate (error_previous_platform, &
                source='     discretization error = 2.664535e-14')
        end select

        !  Test all analysis and synthesis subroutines.
        !  Set scalar field as (x*y*z)**k) restricted to the sphere
        associate (&
            se => original_scalar_function, &
            radial => sphere_type%unit_vectors%radial &
           )
            do k=1, NSYNTHS
                do j=1, NLONS
                    do i=1, NLATS
                        associate (&
                            x => radial(i, j)%x, &
                            y => radial(i, j)%y, &
                            z => radial(i, j)%z &
                           )
                            select case (k)
                                case(1)
                                    se(i, j, k) = exp(x + y + z)
                                case default
                                    se(i, j, k) = (x*y*z)**k
                            end select
                        end associate
                    end do
                end do
            end do
        end associate

        !  Perform analysis then synthesis
        do k = 1, NSYNTHS
            associate (&
                se => original_scalar_function(:, :, k), &
                s => approximate_scalar_function(:, :, k) &
               )

                ! Analyse function into (real) coefficients
                call sphere_type%perform_complex_analysis(se)

                ! Synthesize function from (real) coefficients
                call sphere_type%perform_complex_synthesis(s)
            end associate
        end do

        !  Compute discretization error
        associate (&
            se => original_scalar_function, &
            s => approximate_scalar_function &
           )
            associate (err2 => maxval(abs(s - se)))

                !  Print earlier output from platform with 64-bit floating point
                !    arithmetic followed by the output from this computer
                write( stdout, '(a)') ''
                write( stdout, '(/a/)') '     test analysis and synthesis routines *** TEST RUN *** '
                write( stdout, '(a)') ''
                write( stdout, '(a)') '     grid type = '//sphere_type%grid%grid_type
                write( stdout, '(a)') '     Testing scalar analysis and synthesis'
                write( stdout, '(2(a, i3))') '     nlat = ', NLATS, ' nlon = ', NLONS
                write( stdout, '(a)') '     Previous 64 bit floating point arithmetic result '
                write( stdout, '(a)') error_previous_platform
                write( stdout, '(a)') '     The output from your computer is: '
                write( stdout, '(a, 1pe15.6)') '     discretization error = ', err2
                write( stdout, '(a)') ''
            end associate
        end associate

        !  Release memory
        call sphere_type%destroy()
        deallocate (error_previous_platform)

    end subroutine test_case

end program test_analysis_and_synthesis_routines
