module type_RegularGrid

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_SphericalGrid, only: &
        SphericalGrid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: RegularGrid



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



    ! Declare constructor
    interface RegularGrid
        module procedure regular_grid_constructor
    end interface



contains



    function regular_grid_constructor(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip),         intent (in) :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip),         intent (in) :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        type (RegularGrid)                :: return_value
        !----------------------------------------------------------------------

        call return_value%create(nlat, nlon)

    end function regular_grid_constructor



    subroutine copy_regular_grid(this, object_to_be_copied)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (RegularGrid), intent (out) :: this
        class (RegularGrid), intent (in)  :: object_to_be_copied
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (object_to_be_copied%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (RegularGrid): '&
                //'in assignment (=) '
        end if

        !
        !==> Make copies
        !
        this%initialized = object_to_be_copied%initialized
        this%NUMBER_OF_LONGITUDES = object_to_be_copied%NUMBER_OF_LONGITUDES
        this%NUMBER_OF_LATITUDES = object_to_be_copied%NUMBER_OF_LATITUDES
        this%LONGITUDINAL_MESH = object_to_be_copied%LONGITUDINAL_MESH
        this%LATITUDINAL_MESH = object_to_be_copied%LATITUDINAL_MESH
        this%latitudes = object_to_be_copied%latitudes
        this%longitudes = object_to_be_copied%longitudes
        this%grid_type = object_to_be_copied%grid_type

    end subroutine copy_regular_grid



    subroutine create_regular_grid(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid), intent (in out) :: this
        integer (ip),         intent (in)    :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip),         intent (in)    :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        !
        !==> Set contants
        !
        this%NUMBER_OF_LATITUDES = nlat
        this%NUMBER_OF_LONGITUDES = nlon

        !
        !==> Set the equally-spaced (regular) grid type
        !
        allocate(this%grid_type, source='regular')

        !
        !==> Set longitudinal grid: 0 <= phi <= 2*pi
        !
        call this%get_equally_spaced_longitudes(nlon, this%longitudes)

        !
        !==> Compute equally-spaced latitudes: 0 <= theta <= pi
        !
        call this%get_equally_spaced_latitudes(nlat, this%latitudes)

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_regular_grid



    subroutine destroy_regular_grid(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check initialization flag
        if (.not.this%initialized) then
            return
        end if

        ! Reset constant
        this%LATITUDINAL_MESH = 0.0_wp

        ! Release parent type
        call this%destroy_grid()

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_regular_grid



    subroutine get_equally_spaced_latitudes(this, nlat, theta)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid),    intent (in out) :: this
        integer (ip),           intent (in)     :: nlat  !! number of latitudinal points
        real (wp), allocatable, intent (out)    :: theta(:) !! latitudes: 0 <= theta <= pi
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: k !! counter
        !----------------------------------------------------------------------

        !
        !==> Check validity of input argument
        !
        if (nlat <= 0) then
            error stop 'Object of class (RegularGrid): '&
                //'invalid argument nlat <= 0 in '&
                //'get_equally_spaced_latitudes'
        end if

        !
        !==> Allocate memory
        !
        allocate( theta(nlat) )

        !
        !==> Compute equally spaced latitudinal grid
        !
        associate( dtheta => this%LATITUDINAL_MESH )

            ! Set equally spaced (uniform) mesh size
            dtheta = PI / (nlat-1)

            ! Compute grid
            theta = [ (real(k - 1, kind=wp) * dtheta, k=1, nlat) ]

        end associate

    end subroutine get_equally_spaced_latitudes



    subroutine unformatted_print(this, header)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularGrid), intent (in out) :: this
        character (len=*),   intent (in)     :: header
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
            error stop 'Uninitialized object of class(RegularGrid): '&
                //'in unformatted_print'
        end if

        ! Write latitudes and longitudes
        call this%print_to_unformatted_binary_files(header)

    end subroutine unformatted_print



    subroutine finalize_regular_grid(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (RegularGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_regular_grid



end module type_RegularGrid
