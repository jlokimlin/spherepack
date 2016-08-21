module type_GaussianGrid

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SphericalGrid, only: &
        SphericalGrid

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GaussianGrid


    ! Declare derived data type
    type, extends (SphericalGrid), public ::  GaussianGrid
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
        real (wp), allocatable, public :: gaussian_weights(:)
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
        !----------------------------------------------------------------------
        procedure, public  :: create => create_gaussian_grid
        procedure, public  :: destroy => destroy_gaussian_grid
        procedure, public  :: get_latitudes_and_gaussian_weights
        procedure, public  :: unformatted_print
        generic,   public  :: assignment (=) => copy_gaussian_grid
        procedure, private :: copy_gaussian_grid
        final              :: finalize_gaussian_grid
        !----------------------------------------------------------------------
    end type GaussianGrid



    ! Declare constructor
    interface GaussianGrid
        module procedure gaussian_grid_constructor
    end interface



contains



    function gaussian_grid_constructor(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip), intent (in) :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        type (GaussianGrid)       :: return_value
        !----------------------------------------------------------------------

        call return_value%create(nlat, nlon)

    end function gaussian_grid_constructor



    subroutine copy_gaussian_grid(this, object_to_be_copied)
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        class (GaussianGrid), intent (out) :: this
        class (GaussianGrid), intent (in)  :: object_to_be_copied
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (object_to_be_copied%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (GaussianGrid): '&
                //'in assignment (=) '
        end if

        !
        !==> Make copies
        !
        this%initialized = object_to_be_copied%initialized
        this%NUMBER_OF_LONGITUDES = object_to_be_copied%NUMBER_OF_LONGITUDES
        this%NUMBER_OF_LATITUDES = object_to_be_copied%NUMBER_OF_LATITUDES
        this%LONGITUDINAL_MESH = object_to_be_copied%LONGITUDINAL_MESH
        this%latitudes = object_to_be_copied%latitudes
        this%longitudes = object_to_be_copied%longitudes
        this%gaussian_weights = object_to_be_copied%gaussian_weights
        this%grid_type = object_to_be_copied%grid_type

    end subroutine copy_gaussian_grid



    subroutine create_gaussian_grid(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), target, intent (in out) :: this
        integer (ip),                 intent (in)     :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip),                 intent (in)     :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set contants
        this%NUMBER_OF_LATITUDES = nlat
        this%NUMBER_OF_LONGITUDES = nlon

        ! Set the gaussian grid type
        allocate( this%grid_type, source='gaussian' )

        ! Set longitudinal grid: 0 <= phi <= 2*pi
        call this%get_equally_spaced_longitudes(nlon, this%longitudes)

        ! Set latitudinal grid: 0 <= theta <= pi
        call this%get_latitudes_and_gaussian_weights(nlat, this%latitudes, this%gaussian_weights)

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_gaussian_grid



    subroutine destroy_gaussian_grid(this)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check initialization flag
        if (.not.this%initialized) return

        !
        !==> Release memory
        !
        if (allocated(this%gaussian_weights)) deallocate(this%gaussian_weights)

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
        ! [2]  http://www2.cisl.ucar.edu/resources/legacy/spherepack/documentation#compute_gaussian_latitudes_and_weights.html
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianGrid),   intent (in out) :: this
        integer (ip),           intent (in)     :: nlat     !! number of latitudinal points
        real (wp), allocatable, intent (out)    :: theta(:) !! latitudinal points: 0 <= theta <= pi
        real (wp), allocatable, intent (out)    :: wts(:)   !! gaussian weights
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)  :: error_flag
        integer (ip)  :: dummy_integer !! unused integer variable to maintain backwards compatibility
        real (wp)     :: dummy_real    !! unused double precision variable to maintain backwards compatibility
        !----------------------------------------------------------------------

        ! Check input argument
        if (nlat <= 0) then
            error stop 'Object of class (GaussianGrid): '&
                //'invalid argument nlat <= 0 '&
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
            call compute_gaussian_latitudes_and_weights(nlat, theta, wts, w, lwork, ierror)

        end associate

        !
        !==> Address error flag
        !
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (GaussianGrid) '&
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
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        character (len=*),    intent (in)     :: header
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)  :: file_unit
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.this%initialized) then
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
        ! Dummy arguments
        !----------------------------------------------------------------------
        type (GaussianGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_gaussian_grid



end module type_GaussianGrid
