module type_GaussianGrid

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GaussianGrid

    !----------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------------------------------------
    character (len=250)   :: error_message        !! Probably long enough
    integer (ip)          :: allocate_status      !! To check allocation status
    integer (ip)          :: deallocate_status    !! To check deallocation status
    real (wp), parameter  :: PI = acos( -1.0_wp )
    real (wp), parameter  :: TWO_PI = 2.0_wp * PI
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, public ::  GaussianGrid
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                public :: initialized = .false.
        character (len=3),      public :: grid_type = 'GAU' !! Either gaussian = 'GAU' (default) or equally-spaced = 'REG'
        real (wp),              public :: MESH_PHI = 0.0_wp  !! Uniform mesh in phi
        real (wp),              public :: MESH_THETA = 0.0_wp  !! Only used in 'REG' grid
        real (wp), allocatable, public :: gaussian_weights(:) !! Used for integration, requires 'GAU' for allocation
        real (wp), allocatable, public :: latitudes(:)  !! 0 <= theta <= pi
        real (wp), allocatable, public :: longitudes(:) !! 0 <= phi <= 2*pi
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure,         public  :: create => create_gaussian_grid
        procedure,         public  :: create_equally_spaced_grid
        procedure,         public  :: destroy => destroy_gaussian_grid
        procedure, nopass, public  :: get_gaussian_weights_and_points !! GAQD
        procedure,         private :: get_equally_spaced_longitudes
        procedure,         private :: get_equally_spaced_latitudes
        procedure,         public  :: unformatted_print => print_to_unformatted_binary_files
        final                      :: finalize_gaussian_grid
        !----------------------------------------------------------------------
    end type GaussianGrid


contains


    subroutine create_gaussian_grid( this, nlat, nlon )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        integer (ip),         intent (in)     :: nlat         !! number of latitudinal points
        integer (ip),         intent (in)     :: nlon         !! number of longitudinal points
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set latitudinal grid: 0 <= theta <= pi
        call this%get_gaussian_weights_and_points( &
            nlat, this%latitudes, this%gaussian_weights )

        ! Set longitudinal grid: 0 <= phi <= 2*pi
        call this%get_equally_spaced_longitudes( nlon, this%longitudes )

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_gaussian_grid


    subroutine create_equally_spaced_grid( this, nlat, nlon )
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
        ! Set latitudinal grid: 0 <= theta <= pi
        !

        ! Set the equally-spaced (regular) grid type
        this%grid_type = 'REG'

        ! Compute equally-spaced latitudes: 0 <= theta <= pi
        call this%get_equally_spaced_latitudes( nlat, this%latitudes )

        ! Set longitudinal grid: 0 <= phi <= 2*pi
        call this%get_equally_spaced_longitudes( nlon, this%longitudes )

        ! Set initialization flag
        this%initialized = .true.

    end subroutine create_equally_spaced_grid


    subroutine destroy_gaussian_grid( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check initialization flag
        if ( .not. this%initialized ) return

        ! Reset default grid type
        this%grid_type = 'GAU'

        ! Reset constants
        this%MESH_PHI = 0.0_wp
        this%MESH_THETA = 0.0_wp

        !
        ! Release memory
        !

        ! Check if array is allocated
        if ( allocated( this%gaussian_weights ) ) then

            ! Deallocate array
            deallocate( this%gaussian_weights, stat=deallocate_status )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                error stop 'TYPE (GaussianGrid): '&
                    //'Deallocating GAUSSIAN_WEIGHTS failed in destroy_gaussian_grid'
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%latitudes ) ) then

            ! Deallocate array
            deallocate( this%latitudes, stat=deallocate_status )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                error stop 'TYPE (GaussianGrid):'&
                    //'Deallocating LATITUDES failed in DESTROY_GAUSSIAN_GRID'
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%longitudes ) ) then

            ! Deallocate array
            deallocate( &
                this%longitudes, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
                write( stderr, '(A)' ) 'Deallocating LONGITUDES failed in destroy_gaussian_grid'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_gaussian_grid


    subroutine get_gaussian_weights_and_points( nlat, theta, wts )
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
        !----------------------------------------------------------------------
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
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip),            intent (in)  :: nlat  !! number of latitudinal points
        real (wp), allocatable, intent (out) :: theta(:) !! latitudinal points: 0 <= theta <= pi
        real (wp), allocatable, intent (out) :: wts(:)   !! gaussian weights
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: error_flag
        integer (ip)  :: dummy_integer !! unused integer variable to maintain backwards compatibility
        real (wp)     :: dummy_real    !! unused double precision variable to maintain backwards compatibility
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------
        ! Allocate latitudes (theta)
        !----------------------------------------------------------------------

        ! Check if array is already allocated
        if ( allocated( theta ) ) then

            ! Deallocate if necessary
            deallocate( &
                theta, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
                write( stderr, '(A)' ) 'Deallocating THETA failed in get_GAUSSIAN_WEIGHTS_AND_POINTS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Allocate array
        allocate( &
            theta( 1:nlat ), &
            stat=allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
            write( stderr, '(A)' ) 'Allocating THETA failed in get_GAUSSIAN_WEIGHTS_AND_POINTS'
            write( stderr, '(A)' ) trim( error_message )
        end if

        !----------------------------------------------------------------------
        ! Allocate gaussian weights (wts)
        !----------------------------------------------------------------------

        ! Check if array is already allocated
        if ( allocated( wts ) ) then

            ! Deallocate if necessary
            deallocate( &
                wts, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
                write( stderr, '(A)' ) 'Deallocating WTS failed in get_GAUSSIAN_WEIGHTS_AND_POINTS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Allocate array
        allocate( &
            wts( nlat ), &
            stat=allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
            write( stderr, '(A)' ) 'Allocating WTS failed in get_GAUSSIAN_WEIGHTS_AND_POINTS'
            write( stderr, '(A)' ) trim( error_message )
        end if

        ! Invoke SPHEREPACK
        associate( &
            w => dummy_real, &
            lwork => dummy_integer, &
            ierror => error_flag &
            )
            call gaqd( nlat, theta, wts, w, lwork, ierror )
        end associate

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                write( stderr, '(A)') 'TYPE (GaussianGrid) in GET_GAUSSIAN_WEIGHTS_AND_POINTS'
                write( stderr, '(A, I2)') 'Argument nlat = ', nlat
                write( stderr, '(A)')     'fails to satisfy nlat >= 0'
            case default
                write( stderr, '(A)') 'TYPE (GaussianGrid)'
                write( stderr, '(A)') 'Undetermined error in GET_GAUSSIAN_WEIGHTS_AND_POINTS'
        end select

    end subroutine get_gaussian_weights_and_points


    subroutine get_equally_spaced_latitudes( this, nlat, theta )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid),   intent (in out) :: this
        integer (ip),           intent (in)     :: nlat  !! number of latitudinal points
        real (wp), allocatable, intent (out)    :: theta(:) !! latitudes: 0 <= theta <= pi
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: k     !! counter
        !----------------------------------------------------------------------

        !
        ! Allocate latitudes (theta)
        !

        ! Check if array is already allocated
        if ( allocated( theta ) ) then

            ! Deallocate if necessary
            deallocate( &
                theta, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
                write( stderr, '(A)' ) 'Deallocating THETA failed in get_EQUALLY_SPACED_LATITUDES'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Allocate array
        allocate( &
            theta( nlat ), &
            stat=allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
            write( stderr, '(A)' ) 'Allocating THETA failed in get_EQUALLY_SPACED_LATITUDES'
            write( stderr, '(A)' ) trim( error_message )
        end if

        ! Set equally spaced (uniform) mesh size
        this%MESH_THETA = PI / nlat

        ! Compute latitudinal grid
        associate( dtheta => this%MESH_THETA )
            do  k = 1, nlat
                theta(k) = real(k - 1, kind=wp) * dtheta
            end do
        end associate

    end subroutine get_equally_spaced_latitudes


    subroutine get_equally_spaced_longitudes( this, nlon, phi )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid),   intent (in out) :: this
        integer (ip),           intent (in)     :: nlon !! number of longitudinal points
        real (wp), allocatable, intent (out)     :: phi(:)  !! longitudes: 0 <= phi <= 2*pi
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)   :: l   !! counter
        !----------------------------------------------------------------------

        !
        ! Allocate longitudes
        !

        ! Check if array is allocated
        if ( allocated( phi ) ) then

            ! Deallocate if already allocated
            deallocate( &
                phi, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
                write( stderr, '(A)' ) 'Deallocating PHI failed in get_EQUALLY_SPACED_LONGITUDES'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Allocate array
        allocate( &
            phi( nlon ), &
            stat=allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            write( stderr, '(A)' ) 'TYPE (GaussianGrid)'
            write( stderr, '(A)' ) 'Allocating PHI failed in get_EQUALLY_SPACED_LONGITUDES'
            write( stderr, '(A)' ) trim( error_message )
        end if

        !----------------------------------------------------------------------
        ! Set equally spaced (uniform) mesh size
        !----------------------------------------------------------------------

        this%MESH_PHI = TWO_PI / nlon

        ! Compute longitudinal grid
        associate( dphi => this%MESH_PHI )
            do l = 1, nlon
                phi(l) = real(l - 1, kind=wp) * dphi
            end do
        end associate

    end subroutine get_equally_spaced_longitudes


    subroutine print_to_unformatted_binary_files( this, header )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianGrid), intent (in out) :: this
        character (len=*),    intent (in)     :: header
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: file_unit
        integer (ip)  :: record_length
        !----------------------------------------------------------------------

        ! Write latitudes
        associate( theta => this%latitudes )
            inquire( iolength=record_length ) theta
            open( newunit=file_unit, file=header//'latitudes.dat', status='replace', &
                form='unformatted', access='direct', recl=record_length )
            write( file_unit, rec=1 ) theta
            close( file_unit )
        end associate

        ! Write longitudes
        associate( phi => this%longitudes )
            inquire( iolength=record_length ) phi
            open( newunit=file_unit, file=header//'longitudes.dat', status='replace', &
                form='unformatted', access='direct', recl=record_length )
            write( file_unit, rec=1 ) phi
            close( file_unit )
        end associate

    end subroutine print_to_unformatted_binary_files


    subroutine finalize_gaussian_grid( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (GaussianGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_gaussian_grid


end module type_GaussianGrid
