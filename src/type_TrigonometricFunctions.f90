module type_TrigonometricFunctions

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_SphericalGrid, only: &
        SphericalGrid

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: TrigonometricFunctions

    ! Declare derived data type
    type, public :: TrigonometricFunctions
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                public :: initialized = .false. !! Initialization flag
        integer (ip),           public :: NUMBER_OF_LONGITUDES = 0 !! number of longitudinal points in phi
        integer (ip),           public :: NUMBER_OF_LATITUDES = 0 !! number of latitudinal points in theta
        real (wp), allocatable, public :: sint(:)  !! sin(theta): 0 <= theta <= pi
        real (wp), allocatable, public :: cost(:)  !! cos(theta): 0 <= theta <= pi
        real (wp), allocatable, public :: sinp(:)  !! sin(phi):   0 <=  phi  <= 2*pi
        real (wp), allocatable, public :: cosp(:)  !! cos(phi):   0 <=  phi  <= 2*pi
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public :: create => create_trigonometric_functions
        procedure, public :: destroy => destroy_trigonometric_functions
        final             :: finalize_trigonometric_functions
        !----------------------------------------------------------------------
    end type TrigonometricFunctions


    ! Declare constructor
    interface TrigonometricFunctions
        module procedure trigonometric_functions_constructor
    end interface



contains



    function trigonometric_functions_constructor(grid) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid), intent (in out) :: grid
        type (TrigonometricFunctions)          :: return_value
        !----------------------------------------------------------------------

        call return_value%create(grid)

    end function trigonometric_functions_constructor



    subroutine create_trigonometric_functions(this, grid_type)
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (TrigonometricFunctions), intent (in out) :: this
        class (SphericalGrid),          intent (in)     :: grid_type
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Check if polymorphic argument is usable
        if ( grid_type%initialized .eqv. .false.) then
            error stop 'Object of class (TrigonometricFunctions): '&
                //'initialized polymorphic argument of class (SphericalGrid) '&
                //'in create_trigonometric_functions'
        end if

        ! Allocate memory
        associate( &
            nlat => grid_type%NUMBER_OF_LATITUDES, &
            nlon => grid_type%NUMBER_OF_LONGITUDES &
            )

            ! Set contants
            this%NUMBER_OF_LATITUDES = nlat
            this%NUMBER_OF_LONGITUDES = nlon

            allocate(this%sint(nlat) )
            allocate(this%cost(nlat) )
            allocate(this%sinp(nlon) )
            allocate(this%cosp(nlon) )
        end associate

        ! compute trigonometric functions
        associate( &
            theta => grid_type%latitudes, &
            phi => grid_type%longitudes, &
            sint => this%sint, &
            cost => this%cost, &
            sinp => this%sinp, &
            cosp => this%cosp &
            )

            sint = sin(theta)
            cost = cos(theta)
            sinp = sin(phi)
            cosp = cos(phi)

        end associate

        ! Set flag
        this%initialized = .true.

    end subroutine create_trigonometric_functions



    subroutine destroy_trigonometric_functions(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (TrigonometricFunctions), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if (.not.this%initialized) return

        ! Release memory
        if (allocated(this%sint)) deallocate(this%sint)
        if (allocated(this%cost)) deallocate(this%cost)
        if (allocated(this%sinp)) deallocate(this%sinp)
        if (allocated(this%cosp)) deallocate(this%cosp)

        ! Reset constants
        this%NUMBER_OF_LONGITUDES = 0
        this%NUMBER_OF_LATITUDES = 0

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_trigonometric_functions



    subroutine finalize_trigonometric_functions(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (TrigonometricFunctions), intent (in out)    :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_trigonometric_functions


end module type_TrigonometricFunctions
