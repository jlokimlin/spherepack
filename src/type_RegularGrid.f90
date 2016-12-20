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
    
    type, extends(SphericalGrid), public :: RegularGrid
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
        real(wp), public :: LATITUDINAL_MESH = 0.0_wp
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
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
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip),         intent(in) :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer(ip),         intent(in) :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        type(RegularGrid)                :: return_value
        !----------------------------------------------------------------------

        call return_value%create(nlat, nlon)

    end function regular_grid_constructor



    subroutine copy_regular_grid(self, object_to_be_copied)
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        class(RegularGrid), intent(out) :: self
        class(RegularGrid), intent(in)  :: object_to_be_copied
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (.not.object_to_be_copied%initialized) then
            error stop 'Uninitialized object of class(RegularGrid): '&
                //'in assignment (=) '
        end if

        !
        !  Make copies
        !
        self%initialized = object_to_be_copied%initialized
        self%NUMBER_OF_LONGITUDES = object_to_be_copied%NUMBER_OF_LONGITUDES
        self%NUMBER_OF_LATITUDES = object_to_be_copied%NUMBER_OF_LATITUDES
        self%LONGITUDINAL_MESH = object_to_be_copied%LONGITUDINAL_MESH
        self%LATITUDINAL_MESH = object_to_be_copied%LATITUDINAL_MESH
        self%latitudes = object_to_be_copied%latitudes
        self%longitudes = object_to_be_copied%longitudes
        self%grid_type = object_to_be_copied%grid_type

    end subroutine copy_regular_grid



    subroutine create_regular_grid(self, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(RegularGrid), intent(inout)  :: self
        integer(ip),         intent(in)    :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer(ip),         intent(in)    :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call self%destroy()

        !
        !  Set contants
        !
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon

        !
        !  Set the equally-spaced (regular) grid type
        !
        allocate(self%grid_type, source='regular')

        !
        !  Set longitudinal grid: 0 <= phi <= 2*pi
        !
        call self%get_equally_spaced_longitudes(nlon, self%longitudes)

        !
        !  Compute equally-spaced latitudes: 0 <= theta <= pi
        !
        call self%get_equally_spaced_latitudes(nlat, self%latitudes)

        ! Set initialization flag
        self%initialized = .true.

    end subroutine create_regular_grid



    subroutine destroy_regular_grid(self)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(RegularGrid), intent(inout)  :: self
        !----------------------------------------------------------------------

        ! Check initialization flag
        if (.not.self%initialized) return

        ! Reset constant
        self%LATITUDINAL_MESH = 0.0_wp

        ! Release parent type
        call self%destroy_grid()

        ! Reset flag
        self%initialized = .false.

    end subroutine destroy_regular_grid



    subroutine get_equally_spaced_latitudes(self, nlat, theta)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(RegularGrid),    intent(inout)  :: self
        integer(ip),           intent(in)     :: nlat  !! number of latitudinal points
        real(wp), allocatable, intent(out)    :: theta(:) !! latitudes: 0 <= theta <= pi
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: k !! counter
        !----------------------------------------------------------------------

        !
        !  Check validity of input argument
        !
        if (nlat <= 0) then
            error stop 'Object of class(RegularGrid): '&
                //'invalid argument nlat <= 0 in '&
                //'get_equally_spaced_latitudes'
        end if

        !
        !  Allocate memory
        !
        allocate( theta(nlat) )

        !
        !  Compute equally spaced latitudinal grid
        !
        associate( dtheta => self%LATITUDINAL_MESH )

            ! Set equally spaced (uniform) mesh size
            dtheta = PI / (nlat-1)

            ! Compute grid
            theta = [ (real(k - 1, kind=wp) * dtheta, k=1, nlat) ]

        end associate

    end subroutine get_equally_spaced_latitudes



    subroutine unformatted_print(self, header)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(RegularGrid), intent(inout)  :: self
        character(len=*),   intent(in)     :: header
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(RegularGrid): '&
                //'in unformatted_print'
        end if

        ! Write latitudes and longitudes
        call self%print_to_unformatted_binary_files(header)

    end subroutine unformatted_print



    subroutine finalize_regular_grid(self)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        type(RegularGrid), intent(inout)  :: self
        !----------------------------------------------------------------------

        call self%destroy()

    end subroutine finalize_regular_grid



end module type_RegularGrid
