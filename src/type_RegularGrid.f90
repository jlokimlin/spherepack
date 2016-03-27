module type_RegularGrid

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_Grid, only: &
        SphericalGrid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: RegularGrid

    !----------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------------------------------------
    integer (ip) :: allocate_status !! To check allocation status
    integer (ip) :: deallocate_status !! To check deallocation status
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, extends (SphericalGrid), public :: RegularGrid
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        real (wp), public :: LATITUDINAL_MESH = 0.0_wp
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public  :: create => create_regular_grid
        procedure, public  :: destroy => destroy_regular_grid
        procedure, private :: get_equally_spaced_latitudes
        procedure, public  :: unformatted_print
        final              :: finalize_regular_grid
        !----------------------------------------------------------------------
    end type RegularGrid


contains


    subroutine create_regular_grid( this, nlat, nlon )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid), intent (in out) :: this
        integer (ip),         intent (in)    :: nlat  !! number of latitudinal points
        integer (ip),         intent (in)    :: nlon  !! number of longitudinal points
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set contants
        this%NUMBER_OF_LATITUDES = nlat
        this%NUMBER_OF_LONGITUDES = nlon

        ! Set the equally-spaced (regular) grid type
        allocate( this%grid_type, source='regular' )

        ! Set longitudinal grid: 0 <= phi <= 2*pi
        call this%get_equally_spaced_longitudes( nlon, this%longitudes )

        ! Compute equally-spaced latitudes: 0 <= theta <= pi
        call this%get_equally_spaced_latitudes( nlat, this%latitudes )

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_regular_grid


    subroutine destroy_regular_grid( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check initialization flag
        if ( this%initialized .eqv. .false. ) return

        ! Reset constant
        this%LATITUDINAL_MESH = 0.0_wp

        ! Release parent type
        call this%destroy_grid()

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_regular_grid


    subroutine get_equally_spaced_latitudes( this, nlat, theta )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid),    intent (in out) :: this
        integer (ip),           intent (in)     :: nlat  !! number of latitudinal points
        real (wp), allocatable, intent (out)    :: theta(:) !! latitudes: 0 <= theta <= pi
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)         :: k !! counter
        real (wp), parameter :: PI = acos( -1.0_wp )
        !----------------------------------------------------------------------

        ! Check input argument
        if ( nlat <= 0 ) then
            error stop 'TYPE (Grid): '&
                //'invalid argument NLAT in GET_EQUALLY_SPACED_LATITUDES'
        end if

        ! Release memory
        if (allocated(theta)) deallocate( theta, stat=deallocate_status )
        ! Check allocation status
        if ( deallocate_status /= 0 ) then
            error stop 'TYPE (Grid): '&
                //'Deallocating THETA failed in GET_EQUALLY_SPACED_LATITUDES'
        end if

        ! Allocate memory
        allocate( theta(nlat), stat=allocate_status )
        ! Check allocation status
        if ( allocate_status /= 0 ) then
            error stop 'TYPE (Grid): '&
                //'Allocating PHI failed in GET_EQUALLY_SPACED_LATITUDES'
        end if

        associate( dtheta => this%LATITUDINAL_MESH )
            ! Set equally spaced (uniform) mesh size
            dtheta = PI / (nlat-1)
            ! Compute  equally spaced (uniform) longitudinal grid
            do k = 1, nlat
                theta(k) = -(PI/2) + real(k - 1, kind=wp) * dtheta
            end do
        end associate

    end subroutine get_equally_spaced_latitudes


    subroutine unformatted_print( this, header )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid), intent (in out) :: this
        character (len=*),   intent (in)     :: header
        !----------------------------------------------------------------------

        call this%print_to_unformatted_binary_files( header )

    end subroutine unformatted_print


    subroutine finalize_regular_grid( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (RegularGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_regular_grid


end module type_RegularGrid

