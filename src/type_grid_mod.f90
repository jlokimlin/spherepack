module type_grid_mod

    use, intrinsic :: iso_fortran_env, only: &
        WP     => REAL64, &
        IP     => INT32, &
        stderr => ERROR_UNIT


    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    ! Public derived data type
    public :: grid_t

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    character (len=200)   :: error_message        !! Probably long enough
    integer (IP)          :: allocate_status      !! To check allocation status
    integer (IP)          :: deallocate_status    !! To check deallocation status
    real (WP), parameter  :: PI     = acos( -1.0_WP )
    real (WP), parameter  :: TWO_PI = 2.0_WP * PI
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, public ::  grid_t

        ! All components are public unless stated otherwise
        !---------------------------------------------------------------------------------
        ! Initialization flag
        !---------------------------------------------------------------------------------
        logical                :: initialized  = .false.
        !---------------------------------------------------------------------------------
        ! Grid type: Either gaussian = 'GAU' (default) or equally-spaced = 'REG'
        !---------------------------------------------------------------------------------
        character (len=3)      :: grid_type   = 'GAU'
        !---------------------------------------------------------------------------------
        ! Real constants
        !---------------------------------------------------------------------------------
        real (WP)              :: MESH_PHI    = 0.0_WP  !! Uniform mesh in phi
        real (WP)              :: MESH_THETA  = 0.0_WP  !! Only used in 'REG' grid
        !---------------------------------------------------------------------------------
        ! Allocatable arrays
        !---------------------------------------------------------------------------------
        real (WP), allocatable :: gaussian_weights(:)      !! Used for integration, requires 'GAU' for allocation
        real (WP), allocatable :: latitudes(:)             !! 0 <= theta <= pi
        real (WP), allocatable :: longitudes(:)            !! 0 <= phi <= 2*pi
        !---------------------------------------------------------------------------------

    contains

        ! All methods are private unless stated otherwise
        private

        !---------------------------------------------------------------------------------
        ! Public methods
        !---------------------------------------------------------------------------------
        procedure, non_overridable, public :: Create => Create_gaussian_grid
        procedure, non_overridable, public :: Create_equally_spaced_grid
        procedure, non_overridable, public :: Destroy
        procedure, nopass, public          :: Get_gaussian_weights_and_points !! GAQD
        !---------------------------------------------------------------------------------
        ! Private methods
        !---------------------------------------------------------------------------------
        procedure                          :: Get_equally_spaced_longitudes
        procedure                          :: Get_equally_spaced_latitudes
        !---------------------------------------------------------------------------------
        ! Finalizer
        !---------------------------------------------------------------------------------
        final                              :: Finalize_grid
        !---------------------------------------------------------------------------------

    end type grid_t

contains
    !
    !*****************************************************************************************
    !
    subroutine Create_gaussian_grid( this, nlat, nlon )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), intent (in out)       :: this
        integer (IP),   intent (in)           :: nlat         !! number of latitudinal points
        integer (IP),   intent (in)           :: nlon         !! number of longitudinal points
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        if ( this%initialized ) then

            write( stderr, '(A)' ) 'ERROR: TYPE (grid_t)'
            write( stderr, '(A)' ) 'You must destroy object before calling CREATE_GAUSSIAN_GRID'

        end if

        !--------------------------------------------------------------------------------
        ! Set latitudinal grid: 0 <= theta <= pi
        !--------------------------------------------------------------------------------

        call this%Get_gaussian_weights_and_points( &
            nlat, this%latitudes, this%gaussian_weights )

        !--------------------------------------------------------------------------------
        ! Set longitudinal grid: 0 <= phi <= 2*pi
        !--------------------------------------------------------------------------------

        call this%Get_equally_spaced_longitudes( nlon, this%longitudes )

        !--------------------------------------------------------------------------------
        ! Set initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .true.

    end subroutine Create_gaussian_grid
    !
    !*****************************************************************************************
    !
    subroutine Create_equally_spaced_grid( this, nlat, nlon )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), intent (in out)       :: this
        integer (IP),   intent (in)           :: nlat         !! number of latitudinal points
        integer (IP),   intent (in)           :: nlon         !! number of longitudinal points
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        if ( this%initialized ) then

            write( stderr, '(A)' ) 'TYPE (grid_t)'
            write( stderr, '(A)' ) 'You must destroy object before calling CREATE_EQUALLY_SPACED_GRID'

        end if

        !--------------------------------------------------------------------------------
        ! Set latitudinal grid: 0 <= theta <= pi
        !--------------------------------------------------------------------------------

        ! Set the equally-spaced (regular) grid type
        this%grid_type = 'REG'

        ! Compute equally-spaced latitudes: 0 <= theta <= pi
        call this%Get_equally_spaced_latitudes( nlat, this%latitudes )

        !--------------------------------------------------------------------------------
        ! Set longitudinal grid: 0 <= phi <= 2*pi
        !--------------------------------------------------------------------------------

        call this%Get_equally_spaced_longitudes( nlon, this%longitudes )

        !--------------------------------------------------------------------------------
        ! Set initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .true.

    end subroutine Create_equally_spaced_grid
    !
    !*****************************************************************************************
    !
    subroutine Destroy( this )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        if ( .not. this%initialized ) return

        !--------------------------------------------------------------------------------
        ! Reset default grid type
        !--------------------------------------------------------------------------------

        this%grid_type = 'GAU'

        !--------------------------------------------------------------------------------
        ! Reset constants
        !--------------------------------------------------------------------------------

        this%MESH_PHI   = 0.0_WP
        this%MESH_THETA = 0.0_WP

        !--------------------------------------------------------------------------------
        ! Deallocate arrays
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%gaussian_weights ) ) then

            ! Deallocate array
            deallocate ( &
                this%gaussian_weights, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (grid_t)'
                write( stderr, '(A)' ) 'Deallocating GAUSSIAN_WEIGHTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%latitudes ) ) then

            ! Deallocate array
            deallocate ( &
                this%latitudes, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (grid_t)'
                write( stderr, '(A)' ) 'Deallocating LATITUDES failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%longitudes ) ) then

            ! Deallocate array
            deallocate ( &
                this%longitudes, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (grid_t)'
                write( stderr, '(A)' ) 'Deallocating LONGITUDES failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Reset status
        !--------------------------------------------------------------------------------

        this%initialized = .false.

    end subroutine Destroy
    !
    !*****************************************************************************************
    !
    subroutine Get_gaussian_weights_and_points( nlat, theta, wts )
        !
        !<Purpose:
        !
        ! Computes the nlat-many gaussian (co)latitudes and weights.
        ! the colatitudes are in radians and lie in the interval (0,pi).
        !
        ! References:
        !
        ! [1] Swarztrauber, Paul N.
        !     "On computing the points and weights for Gauss--Legendre quadrature."
        !     SIAM Journal on Scientific Computing 24.3 (2003): 945-954.
        !
        ! [2]  http://www2.cisl.ucar.edu/spherepack/documentation#gaqd.html
        !
        !
        !--------------------------------------------------------------------------------
        ! Documentation: SPHEREPACK 3.2
        !
        !     Gaussian points and weights are computed using the fourier-newton
        !     described in [1]
        !
        !     a test program is located in file tgaqd.f
        !
        !     subroutine gaqd computes the nlat gaussian colatitudes and weights
        !     in double precision. the colatitudes are in radians and lie in the
        !     in the interval (0,pi).
        !
        !     input parameters
        !
        !     nlat    the number of gaussian colatitudes in the interval (0,pi)
        !             (between the two poles).  nlat must be greater than zero.
        !
        !     w       unused double precision variable that permits a simple
        !             exchange with the old routine with the same name in
        !             spherepack.
        !
        !     lwork   unused variable that permits a simple exchange with the
        !             old routine with the same name in SPHEREPACK.
        !
        !     output parameters
        !
        !     theta   a double precision array with length nlat
        !             containing the gaussian colatitudes in
        !             increasing radians on the interval (0,pi).
        !
        !     wts     a double precision array with lenght nlat
        !             containing the gaussian weights.
        !
        !     ierror = 0 no errors
        !            = 1 if nlat <=0
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP),           intent (in)  :: nlat  !! number of latitudinal points
        real (WP), allocatable, intent (out) :: theta(:) !! latitudinal points: 0 <= theta <= pi
        real (WP), allocatable, intent (out) :: wts(:)   !! gaussian weights
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  :: error_flag
        integer (IP)  :: dummy_integer !! unused integer variable to maintain backwards compatibility
        real (WP)     :: dummy_real    !! unused double precision variable to maintain backwards compatibility
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Allocate latitudes (theta)
        !--------------------------------------------------------------------------------

        ! Check if array is already allocated
        if ( allocated( theta ) ) then

            ! Deallocate if necessary
            deallocate ( &
                theta, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (grid_t)'
                write( stderr, '(A)' ) 'Deallocating THETA failed in GET_GAUSSIAN_WEIGHTS_AND_POINTS'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Allocate array
        allocate ( &
            theta( 1:nlat ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (grid_t)'
            write( stderr, '(A)' ) 'Allocating THETA failed in GET_GAUSSIAN_WEIGHTS_AND_POINTS'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Allocate gaussian weights (wts)
        !--------------------------------------------------------------------------------

        ! Check if array is already allocated
        if ( allocated( wts ) ) then

            ! Deallocate if necessary
            deallocate ( &
                wts, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (grid_t)'
                write( stderr, '(A)' ) 'Deallocating WTS failed in GET_GAUSSIAN_WEIGHTS_AND_POINTS'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Allocate array
        allocate ( &
            wts( 1:nlat ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (grid_t)'
            write( stderr, '(A)' ) 'Allocating WTS failed in GET_GAUSSIAN_WEIGHTS_AND_POINTS'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            w      => dummy_real, &
            lwork  => dummy_integer, &
            ierror => error_flag &
            )

            call Gaqd( nlat, theta, wts, w, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0) then

            write( stderr, '(A)') 'ERROR: TYPE (grid_t) in GET_GAUSSIAN_WEIGHTS_AND_POINTS'

        else if (error_flag == 1) then

            write( stderr, '(A, I2)') 'Argument nlat = ', nlat
            write( stderr, '(A)')     'fails to satisfy nlat >= 0'

        else

            write( stderr, '(A)') 'Undetermined error'

        end if

    end subroutine Get_gaussian_weights_and_points
    !
    !*****************************************************************************************
    !
    subroutine Get_equally_spaced_latitudes( this, nlat, theta )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t),         intent (in out) :: this
        integer (IP),           intent (in)     :: nlat  !! number of latitudinal points
        real (WP), allocatable, intent (out)    :: theta(:) !! latitudes: 0 <= theta <= pi
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: k     !! counter
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Allocate latitudes (theta)
        !--------------------------------------------------------------------------------

        ! Check if array is already allocated
        if ( allocated( theta ) ) then

            ! Deallocate if necessary
            deallocate ( &
                theta, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (grid_t)'
                write( stderr, '(A)' ) 'Deallocating THETA failed in GET_EQUALLY_SPACED_LATITUDES'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Allocate array
        allocate ( &
            theta( 1:nlat ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (grid_t)'
            write( stderr, '(A)' ) 'Allocating THETA failed in GET_EQUALLY_SPACED_LATITUDES'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Set equally spaced (uniform) mesh size
        !--------------------------------------------------------------------------------

        this%MESH_THETA = PI / nlat

        !--------------------------------------------------------------------------------
        ! Compute latitudinal grid
        !--------------------------------------------------------------------------------

        associate( Dtheta => this%MESH_THETA )

            do concurrent ( k = 1:nlat )

                theta(k) = real(k - 1, WP) * Dtheta

            end do

        end associate

    end subroutine Get_equally_spaced_latitudes
    !
    !*****************************************************************************************
    !
    subroutine Get_equally_spaced_longitudes( this, nlon, phi )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t),         intent (in out) :: this
        integer (IP),           intent (in)     :: nlon !! number of longitudinal points
        real (WP), allocatable, intent (out)    :: phi(:)  !! longitudes: 0 <= phi <= 2*pi
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                            :: l   !! counter
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Allocate longitudes
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( phi ) ) then

            ! Deallocate if already allocated
            deallocate ( &
                phi, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (grid_t)'
                write( stderr, '(A)' ) 'Deallocating PHI failed in GET_EQUALLY_SPACED_LONGITUDES'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Allocate array
        allocate ( &
            phi( 1:nlon ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (grid_t)'
            write( stderr, '(A)' ) 'Allocating PHI failed in GET_EQUALLY_SPACED_LONGITUDES'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Set equally spaced (uniform) mesh size
        !--------------------------------------------------------------------------------

        this%MESH_PHI = TWO_PI / nlon

        !--------------------------------------------------------------------------------
        ! Compute longitudinal grid
        !--------------------------------------------------------------------------------

        associate( Dphi => this%MESH_PHI )

            do concurrent ( l = 1:nlon )

                phi(l) = real(l - 1, WP) * Dphi

            end do

        end associate

    end subroutine Get_equally_spaced_longitudes
    !
    !*****************************************************************************************
    !
    subroutine Finalize_grid( this )
        !
        ! Purpose:
        !< Finalize object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (grid_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%Destroy()

    end subroutine Finalize_grid
    !
    !*****************************************************************************************
    !
end module type_grid_mod
