module type_GaussianGrid

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_SphericalGrid, only: &
        SphericalGrid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GaussianGrid

    ! Declare derived data type
    type, extends (SphericalGrid), public ::  GaussianGrid
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        real (wp), allocatable, public :: gaussian_weights(:) !! Used for integration, requires 'GAU' for allocation
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public  :: create => create_gaussian_grid
        procedure, public  :: destroy => destroy_gaussian_grid
        procedure, public  :: get_latitudes_and_gaussian_weights
        procedure, public  :: unformatted_print
        final              :: finalize_gaussian_grid
        !----------------------------------------------------------------------
    end type GaussianGrid


contains


    subroutine create_gaussian_grid(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        integer (ip),         intent (in)     :: nlat         !! number of latitudinal points
        integer (ip),         intent (in)     :: nlon         !! number of longitudinal points
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        !
        !==> Set contants
        !
        this%NUMBER_OF_LATITUDES = nlat
        this%NUMBER_OF_LONGITUDES = nlon

        !
        !==> Set the gaussian grid type
        !
        allocate(this%grid_type, source='gaussian')

        !
        !==> Set longitudinal grid: 0 <= phi <= 2*pi
        !
        call this%get_equally_spaced_longitudes(nlon, this%longitudes)

        !
        !==> Set latitudinal grid: 0 <= theta <= pi
        !
        call this%get_latitudes_and_gaussian_weights( &
            nlat, this%latitudes, this%gaussian_weights)

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_gaussian_grid


    subroutine destroy_gaussian_grid(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check initialization flag
        if (this%initialized .eqv. .false.) then
            return
        end if

        !
        !==> Release memory
        !
        if (allocated(this%gaussian_weights)) then
            deallocate(this%gaussian_weights)
        end if

        !
        !==> Release parent type
        !
        call this%destroy_grid()

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_gaussian_grid


    subroutine get_latitudes_and_gaussian_weights(this, nlat, theta, wts)
        !
        !<Purpose:
        !
        ! Computes the nlat-many gaussian (co)latitudes and weights.
        ! the colatitudes are in radians and lie in the interval (0, pi).
        !
        ! References:
        !
        ! [1] Swarztrauber, Paul N.
        !     "On computing the points and weights for Gauss--Legendre quadrature."
        !     SIAM Journal on Scientific Computing 24.3 (2003): 945-954.
        !
        ! [2]  http://www2.cisl.ucar.edu/resources/legacy/spherepack/documentation#gaqd.html
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid),   intent (in out) :: this
        integer (ip),           intent (in)     :: nlat  !! number of latitudinal points
        real (wp), allocatable, intent (out)    :: theta(:) !! latitudinal points: 0 <= theta <= pi
        real (wp), allocatable, intent (out)    :: wts(:)   !! gaussian weights
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: error_flag
        integer (ip)  :: dummy_integer !! unused integer variable to maintain backwards compatibility
        real (wp)     :: dummy_real    !! unused double precision variable to maintain backwards compatibility
        !----------------------------------------------------------------------

        ! Check input argument
        if ( nlat <= 0 ) then
            error stop 'Object of class (GaussianGrid): '&
                //'invalid argument NLAT <= 0 '&
                //'in get_equally_spaced_latitudes'
        end if

        !
        !==> Allocate memory
        !
        allocate( theta(nlat) )
        allocate( wts(nlat) )

        ! Associate various quantities
        associate( &
            w => dummy_real, &
            lwork => dummy_integer, &
            ierror => error_flag &
            )
            !
            !==> Compute gaussian weights and latitudes
            !
            call gaqd(nlat, theta, wts, w, lwork, ierror)
        end associate

        !
        !==> Address error flag
        !
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (GaussianGrid): '&
                    //'fails to satisfy nlat >= 0 '&
                    //'in get_latitudes_and_gaussian_weights'
            case default
                error stop 'Object of class (GaussianGrid): '&
                    //'undetermined error in '&
                    //'get_latitudes_and_gaussian_weights'
        end select

    end subroutine get_latitudes_and_gaussian_weights


    subroutine unformatted_print(this, header)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        character (len=*),    intent (in)     :: header
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: file_unit
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (GaussianGrid): '&
                //'in unformatted_print'
        end if

        ! Write latitudes and longitudes
        call this%print_to_unformatted_binary_files(header)

        ! Write gaussian weights
        associate( wts => this%gaussian_weights )

            open( newunit=file_unit, file=header//'gaussian_weights.dat', &
                status='replace', form='unformatted', &
                action='write', access='stream' )
            write( file_unit ) wts
            close( file_unit )

        end associate

    end subroutine unformatted_print


    subroutine finalize_gaussian_grid(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (GaussianGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_gaussian_grid


end module type_GaussianGrid
