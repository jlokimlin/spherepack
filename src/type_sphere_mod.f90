!*****************************************************************************************
!
!< Author:
! Jon Lo Kim Lin
!
!< Purpose:
!
! A modern Fortran (2008+) object-oriented spherepack wrapper for NCAR's SPHEREPACK 3.2
!
!*****************************************************************************************
!
module type_sphere_mod

    use, intrinsic :: iso_fortran_env, only: &
        WP     => REAL64, &
        IP     => INT32, &
        stderr => ERROR_UNIT

    use type_workspace_mod, only: &
        workspace_t

    use type_grid_mod, only: &
        grid_t

    use type_vector_mod, only: &
        vector_t, &
        assignment(=), &
        operator(*)
    
    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sphere_t

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    character (len=200)   :: error_message     !! Probably long enough
    integer (IP)          :: allocate_status   !! To check allocation status
    integer (IP)          :: deallocate_status !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, public :: sphere_t

        ! All components are public unless stated otherwise

        !---------------------------------------------------------------------------------
        ! Initialization flag
        !---------------------------------------------------------------------------------
        logical                      :: initialized = .false. !! Instantiation status
        !---------------------------------------------------------------------------------
        ! Integer constants
        !---------------------------------------------------------------------------------
        integer (IP)                 :: NLON                = 0   !! number of longitudinal points
        integer (IP)                 :: NLAT                = 0   !! number of latitudinal points
        integer (IP)                 :: NTRUNC              = 0   !! triangular truncation limit
        integer (IP)                 :: SCALAR_SYMMETRIES   = 0   !! symmetries about the equator for scalar calculations
        integer (IP)                 :: VECTOR_SYMMETRIES   = 0   !! symmetries about the equator for vector calculations
        integer (IP)                 :: NUMBER_OF_SYNTHESES = 0   !!
        !---------------------------------------------------------------------------------
        ! Complex spectral coefficients
        !---------------------------------------------------------------------------------
        complex (WP), allocatable    :: complex_spectral_coefficients(:)
        !---------------------------------------------------------------------------------
        ! Derived data type
        !---------------------------------------------------------------------------------
        type (workspace_t), private  :: workspace             !! Contains the various workspace arrays to invoke SPHERPACK 3.2
        type (grid_t)                :: grid                  !! Spherical grid
        !---------------------------------------------------------------------------------
        ! Commonly used trigonometric functions
        !---------------------------------------------------------------------------------
        real (WP), allocatable       :: sint(:)                 !! sin(theta): 0 <= theta <= pi
        real (WP), allocatable       :: cost(:)                 !! cos(theta): 0 <= theta <= pi
        real (WP), allocatable       :: sinp(:)                 !! sin(phi):   0 <=  phi  <= 2*pi
        real (WP), allocatable       :: cosp(:)                 !! cos(phi):   0 <=  phi  <= 2*pi
        !---------------------------------------------------------------------------------
        ! The spherical unit vectors
        !---------------------------------------------------------------------------------
        real (WP), allocatable       :: radial_unit_vector(:, :, :)
        real (WP), allocatable       :: polar_unit_vector(:, :, :)
        real (WP), allocatable       :: azimuthal_unit_vector(:, :, :)
        !---------------------------------------------------------------------------------

    contains
        
        ! All method are private unless stated otherwise
        private

        !---------------------------------------------------------------------------------
        ! Public SPHEREPACK 3.2 methods
        !---------------------------------------------------------------------------------
        procedure, public                    :: Get_colatitude_derivative !! Vtsgs
        procedure, public                    :: Get_Gradient !! Gradgs
        procedure, public                    :: Invert_gradient !!  Igradgs
        procedure, public                    :: Get_Divergence !! Divgs
        procedure, public                    :: Invert_divergence !!Idivgs
        procedure, public                    :: Get_Vorticity !! Vrtgs
        procedure, public                    :: Invert_vorticity !! Ivrtgs
        procedure, public                    :: Invert_divergence_and_vorticity !! Idvtgs
        procedure, public                    :: Get_Scalar_laplacian !! Slapgs
        procedure, public                    :: Invert_helmholtz !! Islapgs
        procedure, public                    :: Get_Vector_laplacian !! Vlapgs
        procedure, public                    :: Invert_vector_laplacian !! Ivlapgs
        procedure, public                    :: Get_Stream_function_and_velocity_potential
        procedure, public                    :: Invert_stream_function_and_velocity_potential
        procedure, public                    :: Perform_grid_transfers
        procedure, public                    :: Perform_Geo_math_coordinate_transfers
        procedure, public                    :: Perform_scalar_analysis
        procedure, public                    :: Perform_scalar_synthesis
        procedure, public                    :: Perform_scalar_projection !! Shpg
        procedure, public                    :: Perform_vector_analysis
        procedure, public                    :: Perform_vector_synthesis
        procedure, public                    :: Get_Legendre_functions
        !procedure, public                    :: Get_Icosahedral_geodesic
        !procedure, public                    :: Get_Multiple_ffts
        procedure, nopass, public            :: Get_gaussian_weights_and_points !! Gaqd
        !procedure, public                    :: Get_Three_dimensional_sphere_graphics
        !---------------------------------------------------------------------------------
        ! Public complex methods
        !---------------------------------------------------------------------------------
        procedure, public                    :: Perform_complex_analysis
        procedure, public                    :: Perform_complex_synthesis
        !---------------------------------------------------------------------------------
        ! Additional public methods
        !---------------------------------------------------------------------------------
        procedure, non_overridable, public   :: Create
        procedure, non_overridable, public   :: Destroy
        procedure, non_overridable, public   :: Get_index
        procedure, non_overridable, public   :: Get_coefficient
        procedure, public                    :: Compute_surface_integral
        procedure, public                    :: Get_rotation_operator
        procedure, public                    :: Synthesize_from_spec
        !---------------------------------------------------------------------------------
        ! Private methods
        !---------------------------------------------------------------------------------
        procedure, non_overridable           :: Assert_initialized
        procedure                            :: Set_trigonometric_functions
        procedure                            :: Set_spherical_unit_vectors
        procedure                            :: Get_spherical_angle_components
        procedure                            :: Set_scalar_symmetries
        procedure                            :: Set_vector_symmetries
        final                                :: Finalize_sphere
        !---------------------------------------------------------------------------------

    end type sphere_t

contains
    !
    !*****************************************************************************************
    !
    subroutine Create( this, nlat, nlon, isym, itype )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        integer (IP),     intent (in)            :: nlat
        integer (IP),     intent (in)            :: nlon
        integer (IP),     intent (in), optional  :: isym      !! Either 0, 1, or 2
        integer (IP),     intent (in), optional  :: itype     !! Either 0, 1, 2, 3, ..., 8
        !character (*), intent (in), optional :: grid_type !! Either '(GAU)' or '(REG)'
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        if ( this%initialized ) then

            write( stderr, '(A)' ) 'ERROR: TYPE (sphere_t)'
            write( stderr, '(A)' ) 'You must destroy object before calling CREATE'

        end if

        !--------------------------------------------------------------------------------
        ! Set constants
        !--------------------------------------------------------------------------------

        this%NLAT                = nlat
        this%NLON                = nlon
        this%NTRUNC              = nlat - 1 !! Set triangular truncation
        this%NUMBER_OF_SYNTHESES = 1

        ! Set scalar symmetries
        if ( present( isym ) ) then

            call this%Set_scalar_symmetries( isym )

        end if

        ! Set vector symmetries
        if (present (itype) ) then

            call this%Set_vector_symmetries( itype )

        end if

        !--------------------------------------------------------------------------------
        ! Allocate array
        !--------------------------------------------------------------------------------

        associate( size_spec => nlat * (nlat + 1)/2 )

            ! Allocate pointer for complex spectral coefficients
            allocate ( &
                this%complex_spectral_coefficients( 1:size_spec ), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocate status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Allocating COMPLEX_SPECTRAL_COEFFICIENTS failed in CREATE'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! Create derived data types
        !--------------------------------------------------------------------------------

        call this%grid%Create( nlat, nlon )
        call this%workspace%Create( nlat, nlon )


        !--------------------------------------------------------------------------------
        ! Set frequently used trigonometric functions
        !--------------------------------------------------------------------------------

        call this%Set_trigonometric_functions( &
            this%grid%latitudes, &
            this%grid%longitudes )

        !--------------------------------------------------------------------------------
        ! Set spherical unit vectors to compute polar and azimuthal components for vector functions
        !--------------------------------------------------------------------------------

        call this%Set_spherical_unit_vectors( &
            this%sint, &
            this%cost, &
            this%sinp, &
            this%cosp )

        !--------------------------------------------------------------------------------
        ! Set initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .true.
        
    end subroutine Create
    !
    !*****************************************************************************************
    !
    subroutine Destroy( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check status
        !--------------------------------------------------------------------------------

        if ( .not. this%initialized ) return

        !--------------------------------------------------------------------------------
        ! Reset constants
        !--------------------------------------------------------------------------------

        this%NLON                = 0
        this%NLAT                = 0
        this%NTRUNC              = 0
        this%SCALAR_SYMMETRIES   = 0
        this%VECTOR_SYMMETRIES   = 0
        this%NUMBER_OF_SYNTHESES = 0

        !--------------------------------------------------------------------------------
        ! Deallocate complex spectral coefficients
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%complex_spectral_coefficients ) ) then

            ! Deallocate array
            deallocate( &
                this%complex_spectral_coefficients, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating COMPLEX_SPECTRAL_COEFFICIENTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end if

        !--------------------------------------------------------------------------------
        ! Destroy derived data types
        !--------------------------------------------------------------------------------

        call this%grid%Destroy()
        call this%workspace%Destroy()

        !--------------------------------------------------------------------------------
        ! Deallocate trigonometric functions
        !--------------------------------------------------------------------------------

         ! Check if array is allocated
        if ( allocated( this%sint ) ) then

            ! Deallocate array
            deallocate( &
                this%sint, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating SINT failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%cost ) ) then

            ! Deallocate array
            deallocate( &
                this%cost, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating COST failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%sinp ) ) then

            ! Deallocate array
            deallocate( &
                this%sinp, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating SINP failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end if

        ! Check if array is allocated
        if ( allocated( this%cosp ) ) then

            ! Deallocate array
            deallocate( &
                this%cosp, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating COSP failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end if

        !--------------------------------------------------------------------------------
        ! Deallocate spherical unit vectors
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%radial_unit_vector ) ) then

            ! Deallocate array
            deallocate( &
                this%radial_unit_vector, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating RADIAL_UNIT_VECTOR failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%polar_unit_vector ) ) then

            ! Deallocate array
            deallocate( &
                this%polar_unit_vector, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating POLAR_UNIT_VECTOR failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if
        
        ! Check if array is allocated
        if ( allocated( this%azimuthal_unit_vector ) ) then

            ! Deallocate array
            deallocate( &
                this%azimuthal_unit_vector, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Deallocating AZIMUTHAL_UNIT_VECTOR failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Reset initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .false.

    end subroutine Destroy
    !
    !*****************************************************************************************
    !
    subroutine Assert_initialized( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. this%initialized ) then
            error stop 'ERROR: You must instantiate type(sphere_t) before calling methods'
        end if

    end subroutine Assert_initialized
    !
    !*****************************************************************************************
    !
    function Get_index( this, n, m ) result( return_value )
        !
        !< Purpose:
        ! The spectral data is assumed to be in a complex array of dimension
        ! (MTRUNC+1)*(MTRUNC+2)/2. MTRUNC is the triangular truncation limit
        ! (MTRUNC = 42 for T42). MTRUNC must be <= nlat-1.
        !
        ! The coefficients are ordered so that
        !
        ! first (nm=1)  is m=0, n=0, second (nm=2) is m=0, n=1,
        ! nm=MTRUNC is m=0, n=MTRUNC, nm=MTRUNC+1 is m=1, n=1, etc.
        !
        ! In other words,
        !
        ! 00, 01, 02, 03, 04.........0MTRUNC
        !     11, 12, 13, 14.........1MTRUNC
        !         22, 23, 24.........2MTRUNC
        !             33, 34.........3MTRUNC
        !                 44.........4MTRUNC
        !                    .
        !                      .
        !                        .etc...
        !
        ! In modern Fortran syntax, values of m (degree) and n (order)
        ! as a function of the index nm are:
        !
        ! integer (IP), dimension ((MTRUNC+1)*(MTRUNC+2)/2) :: indxm, indxn
        ! indxm = [((m, n=m, MTRUNC), m=0, MTRUNC)]
        ! indxn = [((n, n=m, MTRUNC), m=0, MTRUNC)]
        !
        ! Conversely, the index nm as a function of m and n is:
        ! nm = sum([(i, i=MTRUNC+1, MTRUNC-m+2, -1)])+n-m+1
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        integer (IP),     intent (in)      :: n
        integer (IP),     intent (in)      :: m
        integer (IP)                       :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: i !! Counter
        !--------------------------------------------------------------------------------

        associate( ntrunc => this%NTRUNC )

            if ( m <= n .and. max( n, m ) <= ntrunc ) then

                return_value = &
                    sum ( [ (i, i = ntrunc+1, ntrunc - m + 2, - 1) ] ) + n - m + 1
            else

                return_value = -1

            end if

        end associate

    end function Get_index
    !
    !*****************************************************************************************
    !
    function Get_coefficient( this, n, m ) result ( return_value )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        integer (IP),     intent (in)        :: n
        integer (IP),     intent (in)        :: m
        complex (WP)                         :: return_value
        !--------------------------------------------------------------------------------

        associate( &
            ntrunc   => this%NTRUNC, &
            nm       => this%Get_index( n, m ), &
            nm_conjg => this%Get_index(n, -m ), &
            psi      => this%complex_spectral_coefficients &
            )

            if ( m < 0 .and. nm_conjg > 0 ) then

                return_value = ( (-1.0_WP)**(-m) ) * conjg( psi(nm_conjg) )

            else if ( nm > 0 ) then

                return_value = psi(nm)

            else

                return_value = 0.0_WP

            end if

        end associate

    end function Get_coefficient
    !
    !*****************************************************************************************
    !
    subroutine Set_scalar_symmetries( this, isym )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        integer (IP),     intent (in)     :: isym
        !--------------------------------------------------------------------------------

        if ( isym == 2 ) then

            this%SCALAR_SYMMETRIES = isym

        else if ( isym == 1) then

            this%SCALAR_SYMMETRIES = isym

        else if ( isym == 0 ) then

            this%SCALAR_SYMMETRIES = isym

        else

            ! Handle invalid symmetry arguments
            write( stderr, '(A)' )     'TYPE (sphere_t)'
            write( stderr, '(A, I2)' ) 'Optional argument isym = ', isym
            write( stderr, '(A)' )     'in SET_SCALAR_SYMMETRIES'
            write( stderr, '(A)' )     'must be either 0, 1, or 2 (default isym = 0)'

        end if

    end subroutine Set_scalar_symmetries
        !
    !*****************************************************************************************
    !
    subroutine Set_vector_symmetries( this, itype )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        integer (IP), intent (in)         :: itype
        !--------------------------------------------------------------------------------

        if ( itype == 8 ) then

            this%VECTOR_SYMMETRIES = itype
            return

        else if ( itype == 7) then

            this%VECTOR_SYMMETRIES = itype
            return

        else if ( itype == 6) then

            this%VECTOR_SYMMETRIES = itype
            return

        else if ( itype == 5) then

            this%VECTOR_SYMMETRIES = itype
            return

        else if ( itype == 4) then

            this%VECTOR_SYMMETRIES = itype
            return

        else if ( itype == 3) then

            this%VECTOR_SYMMETRIES = itype
            return

        else if ( itype == 2) then

            this%VECTOR_SYMMETRIES = itype
            return

        else if ( itype == 1) then

            this%VECTOR_SYMMETRIES = itype
            return


        else if ( itype == 0 ) then

            this%VECTOR_SYMMETRIES = itype
            return

        else

            ! Handle invalid symmetry arguments
            write( stderr, '(A)' )     'TYPE (sphere_t)'
            write( stderr, '(A, I2)' ) 'Optional argument itype = ', itype
            write( stderr, '(A)' )     'in SET_VECTOR_SYMMETRIES'
            write( stderr, '(A)' )     'must be either 0, 1, 2, ..., 8 (default itype = 0)'

        end if

    end subroutine Set_vector_symmetries
    !
    !*****************************************************************************************
    !
    subroutine Set_trigonometric_functions( this, theta, phi )
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        real (WP),        intent (in)     :: theta(:)
        real (WP),        intent (in)     :: phi(:)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if latitudes are allocated
        !--------------------------------------------------------------------------------

        if ( .not. allocated( this%grid%latitudes ) ) then

            write( stderr, '(A)' ) 'TYPE (sphere_t)'
            write( stderr, '(A)' ) 'in SET_TRIGONOMETRIC_FUNCTIONS'
            write( stderr, '(A)' ) 'You must allocate LATITUDES '
            write( stderr, '(A)' ) 'before calling SET_TRIGONOMETRIC_FUNCTIONS'

        end if

        !--------------------------------------------------------------------------------
        ! Check if longitudes are allocated
        !--------------------------------------------------------------------------------

        if ( .not. allocated( this%grid%longitudes ) ) then

            write( stderr, '(A)' ) 'TYPE (sphere_t)'
            write( stderr, '(A)' ) 'in SET_TRIGONOMETRIC_FUNCTIONS'
            write( stderr, '(A)' ) 'You must allocate LONGITUDES '
            write( stderr, '(A)' ) 'before calling SET_TRIGONOMETRIC_FUNCTIONS'

        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        associate( &
            nlat => size( theta ), &
            nlon => size( phi ) &
            )

            allocate ( &
                this%sint( 1:nlat ), &
                this%cost( 1:nlat ), &
                this%sinp( 1:nlon ), &
                this%cosp( 1:nlon ), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)') 'TYPE (sphere_t)'
                write( stderr, '(A)') 'Allocation failed in SET_TRIGONOMETRIC_FUNCTIONS'
                write( stderr, '(A)') trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! Compute trigonometric functions
        !--------------------------------------------------------------------------------

        this%sint = sin( theta )
        this%cost = cos( theta )
        this%sinp = sin( phi )
        this%cosp = cos( phi )

    end subroutine Set_trigonometric_functions
    !
    !*****************************************************************************************
    !
    subroutine Set_spherical_unit_vectors( this, sint, cost, sinp, cosp )
        !
        !< Purpose:
        ! Sets the spherical unit vectors
        !
        ! Remark:
        ! The "grid" component of sphere must be
        ! initialized before calling this subroutine
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)   :: this
        real (WP),        intent (in)       :: sint(:)
        real (WP),        intent (in)       :: cost(:)
        real (WP),        intent (in)       :: sinp(:)
        real (WP),        intent (in)       :: cosp(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  k, l !! Counters
        !--------------------------------------------------------------------------------

        dimensions: associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Allocate arrays
            allocate ( &
                this%radial_unit_vector( 1:3, 1:nlat, 1:nlon ), &
                this%polar_unit_vector( 1:3, 1:nlat, 1:nlon ), &
                this%azimuthal_unit_vector( 1:3, 1:nlat, 1:nlon ), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocate status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)') 'TYPE (sphere_t)'
                write( stderr, '(A)') 'Allocation failed in SET_SPHERICAL_UNIT_VECTORS'
                write( stderr, '(A)') trim( error_message )

            end if

            unit_vectors: associate( &
                r     => this%radial_unit_vector, &
                theta => this%polar_unit_vector, &
                phi   => this%azimuthal_unit_vector &
                )

                ! Compute spherical unit vectors
                do l = 1, nlon
                    do k = 1, nlat

                        ! set radial unit vector
                        r(:, k, l) = &
                            [ sint(k) * cosp(l), &
                            sint(k) * sinp(l), &
                            cost(k) ]

                        ! set polar unit vector
                        theta(:, k, l) = &
                            [ cost(k) * cosp(l), &
                            cost(k) * sinp(l), &
                            -sint(k) ]

                        ! set azimuthal unit vector
                        phi(:, k, l) = &
                            [ -sinp(l), &
                            cosp(l), &
                            0.0_WP ]
                    end do
                end do

            end associate unit_vectors

        end associate dimensions

    end subroutine Set_spherical_unit_vectors
    !
    !*****************************************************************************************
    !
    subroutine Perform_complex_analysis( this, scalar_function )
        !
        !<Purpose:
        !
        ! Converts gridded input array (scalar_function) to (complex)
        ! spherical harmonic coefficients (psi).
        !
        ! The spectral data is assumed to be in a complex array of dimension
        ! (mtrunc+1)*(mtrunc+2)/2, whre mtrunc is the triangular truncation limit,
        ! for instance, mtrunc = 42 for T42.
        ! mtrunc must be <= nlat-1.
        ! Coefficients are ordered so that first (nm=1) is m=0, n=0, second is m=0, n=1,
        ! nm=mtrunc is m=0, n=mtrunc, nm=mtrunc+1 is m=1, n=1, etc.
        !
        ! In modern Fortran syntax, values of m (degree) and n (order) as a function
        ! of the index nm are:
        !
        ! integer (IP), dimension ((mtrunc+1)*(mtrunc+2)/2) :: indxm, indxn
        ! indxm = [((m, n=m, mtrunc), m=0, mtrunc)]
        ! indxn = [((n, n=m, mtrunc), m=0, mtrunc)]
        !
        ! Conversely, the index nm as a function of m and n is:
        ! nm = sum([(i, i=mtrunc+1, mtrunc-m+2, -1)])+n-m+1
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (in)      :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: m,  n  !< counters
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()
        
        !--------------------------------------------------------------------------------
        ! Compute the (real) spherical harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%Perform_scalar_analysis( scalar_function )

        !--------------------------------------------------------------------------------
        ! Set complex spherical harmonic coefficients
        !--------------------------------------------------------------------------------

        associate( &
            ntrunc => this%NTRUNC, & ! set the triangular truncation limit
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients, &
            psi    => this%complex_spectral_coefficients &
            )

            ! Fill complex array
            psi = cmplx( &
                0.5_WP * [((a(m + 1, n + 1), n = m, ntrunc), m = 0, ntrunc)], &
                0.5_WP * [((b(m + 1, n + 1), n = m, ntrunc), m = 0, ntrunc)], &
                WP )
 
        end associate

    end subroutine Perform_complex_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_complex_synthesis( this, scalar_function )
        !
        !< Purpose:
        ! converts gridded input array (datagrid) to (complex) spherical harmonic coefficients
        ! (dataspec).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (out)     :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: m, n  !! Counters
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Convert complex spherical harmonic coefficients to real version
        !--------------------------------------------------------------------------------
        associate( &
            ntrunc => this%NTRUNC, & ! set the triangular truncation limit
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients, &
            psi    => this%complex_spectral_coefficients &
            )
 
            ! Fill real arrays with contents of spec
            do m = 0, ntrunc
                do n = m, ntrunc
                
                    ! set the spectral index
                    associate(nm => this%Get_index( n, m ))
                
                        ! set the real component
                        a( m + 1, n + 1 ) = 2.0_WP * real( psi(nm) )
                
                        ! set the imaginary component
                        b( m + 1, n + 1 ) = 2.0_WP * aimag( psi(nm) )

                    end associate
                end do
            end do
        
        end associate

        !--------------------------------------------------------------------------------
        ! Synthesise the scalar function from the (real) harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%Perform_scalar_synthesis( scalar_function )
 
    end subroutine Perform_complex_synthesis
    !
    !*****************************************************************************************
    !
    subroutine Synthesize_from_spec( this, spec, scalar_function )
        !
        !< Purpose:
        !
        ! Used mainly for testing
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        complex (WP),     intent (in)      :: spec(:)
        real (WP),        intent (out)     :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                       :: m, n  !! Counters
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Convert complex coefficients to real version
        !--------------------------------------------------------------------------------

        associate( &
            ntrunc => this%NTRUNC, & ! set the triangular truncation limit
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients &
            )
 
            ! fill real arrays with contents of spec
            do m = 0, ntrunc
                do n = m, ntrunc
                
                    ! set the spectral index
                    spec_index: associate( nm => this%Get_index( n, m ) )
                
                        ! set the real component
                        a(m + 1, n + 1) = 2.0_WP * real( spec(nm) )
                
                        ! set the imaginary component
                        b(m + 1, n + 1) = 2.0_WP * aimag( spec(nm) )

                    end associate spec_index

                end do
            end do
        
        end associate

        !--------------------------------------------------------------------------------
        ! Synthesise the scalar function from the (real) harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%Perform_scalar_synthesis( scalar_function )
 
    end subroutine Synthesize_from_spec
    !
    !*****************************************************************************************
    !
    subroutine Get_spherical_angle_components( this, &
        vector_function, polar_component, azimuthal_component )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (in)      :: vector_function(:, :, :)
        real (WP),        intent (out)     :: polar_component(:, :)
        real (WP),        intent (out)     :: azimuthal_component(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: k,  l        !! Counters
        type (vector_t) :: theta        !! Polar unit vector
        type (vector_t) :: phi          !! Azimuthal unit vector
        type (vector_t) :: vector_field !! To convert array to vector
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Initialize arrays
        !--------------------------------------------------------------------------------
        
        polar_component     = 0.0_WP
        azimuthal_component = 0.0_WP
        
        !--------------------------------------------------------------------------------
        ! Calculate the spherical angle components
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )


            do l = 1, nlon
                do k = 1, nlat

                    ! Convert arrays to vectors
                    vector_field = vector_function(:, k, l)
                    theta = this%polar_unit_vector(:, k, l)
                    phi = this%azimuthal_unit_vector(:, k, l)

                    ! set the theta component
                    polar_component(k, l) = theta.dot.vector_field
               
                    ! set the azimuthal_component
                    azimuthal_component(k, l) = phi.dot.vector_field

                end do
            end do

        end associate

    end subroutine Get_spherical_angle_components
    !
    !*****************************************************************************************
    !
    subroutine Get_rotation_operator( this, scalar_function, rotation_operator)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (in)      :: scalar_function(:, :)
        real (WP),        intent (out)     :: rotation_operator(:, :, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)           :: k, l   !! Counters
        type (vector_t)        :: theta  !! Polar unit vector
        type (vector_t)        :: phi    !! Azimuthal unit vector
        real (WP), allocatable :: polar_gradient_component(:, :)
        real (WP), allocatable :: azimuthal_gradient_component(:, :)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Allocate arrays
            allocate ( &
                polar_gradient_component(     1:nlat, 1:nlon ), &
                azimuthal_gradient_component( 1:nlat, 1:nlon ), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Allocation failed in GET_ROTATION_OPERATOR'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! Calculate the spherical surface gradient components
        !--------------------------------------------------------------------------------

        associate( &
            f          => scalar_function, &
            grad_theta => polar_gradient_component, &
            grad_phi   => azimuthal_gradient_component &
            )

            call this%Get_gradient( f, grad_theta, grad_phi )

        end associate

        !--------------------------------------------------------------------------------
        ! Calculate the rotation operator applied to a scalar function
        !--------------------------------------------------------------------------------

        associate( &
            nlat       => this%NLAT, &
            nlon       => this%NLON, &
            R          => rotation_operator &
            )

            ! Initialize array
            R = 0.0_WP

            do l = 1, nlon
                do k = 1, nlat

                    ! Convert arrays to vectors
                    theta = this%polar_unit_vector(:, k, l)
                    phi = this%azimuthal_unit_vector(:, k, l)

                    associate( &
                        grad_theta => polar_gradient_component(k, l), &
                        grad_phi   => azimuthal_gradient_component(k, l) &
                        )

                        R(:, k, l) = phi * grad_theta - theta * grad_phi

                    end associate
                end do
            end do
        end associate

        !--------------------------------------------------------------------------------
        ! Deallocate arrays
        !--------------------------------------------------------------------------------

        deallocate ( &
            polar_gradient_component, &
            azimuthal_gradient_component, &
            stat   = deallocate_status, &
            errmsg = error_message )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (sphere_t)'
            write( stderr, '(A)' ) 'Deallocation failed in GET_ROTATION_OPERATOR'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end subroutine Get_rotation_operator
    !
    !*****************************************************************************************
    !
    function Compute_surface_integral( this, scalar_function ) result( return_value )
        !
        !< Purpose:
        !
        ! Computes the (scalar) surface integral on the sphere (S^2):
        !
        ! * Trapezoidal rule    in phi:   0 <=  phi  <= 2*pi
        ! * Gaussian quadrature in theta: 0 <= theta <= pi
        !
        !   \int_{S^2} f( theta, phi ) dS
        !
        !   where
        !
        !   dS = sin(theta) dtheta dphi
        !
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        real (WP),        intent (in)     :: scalar_function(:, :)
        real (WP)                         :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)           :: k            !! counter
        real (WP), allocatable :: summation(:)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate array
        !--------------------------------------------------------------------------------

        associate( nlat => this%NLAT )

            ! Allocate array
            allocate ( &
                summation( 1:nlat ), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Allocation failed in COMPUTE_SURFACE_INTEGRAL'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! Compute integral
        !--------------------------------------------------------------------------------

        ! Initialize array
        summation = 0.0_WP

        ! Compute the integrant
        associate( &
            nlat => this%NLAT, &
            Dphi => this%grid%MESH_PHI, &
            wts  => this%grid%gaussian_weights, &
            f    => scalar_function &
            )

            ! Apply trapezoidal rule
            do k = 1, nlat

                summation(k) = sum( f(k, :) ) * Dphi

            end do

            ! Apply gaussian quadrature
            summation = summation * wts

        end associate

        ! Set integral \int_{S^2} f( theta, phi ) dS
        return_value = sum( summation )

        !--------------------------------------------------------------------------------
        ! Deallocate array
        !--------------------------------------------------------------------------------

        deallocate ( &
            summation, &
            stat   = deallocate_status, &
            errmsg = error_message )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (sphere_t)'
            write( stderr, '(A)' ) 'Deallocation failed in COMPUTE_SURFACE_INTEGRAL'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end function Compute_surface_integral
    !
    !*****************************************************************************************
    !
    subroutine Compute_first_moment( this, scalar_function, first_moment )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (in)      :: scalar_function(:, :)
        type (vector_t),  intent (out)     :: first_moment
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)           :: k,  l             !! Counters
        type (vector_t)        :: u                 !! radial unit vector
        real (WP), allocatable :: integrant(:, :, :)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate array
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Allocate arrays
            allocate ( &
                integrant( 1:nlat, 1:nlon, 1:3 ), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Allocation failed in COMPUTE_FIRST_MOMENT'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! Compute integrant
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Initialize array
            integrant = 0.0_WP

            ! Compute integrant
            do l = 1, nlon
                do k = 1, nlat

                    ! Convert array to vector
                    u = this%radial_unit_vector(:, k, l)

                    associate( f => scalar_function(k, l) )

                        integrant(k, l, 1) = u%x * f
                        integrant(k, l, 2) = u%y * f
                        integrant(k, l, 3) = u%z * f

                    end associate
                end do
            end do

        end associate

        !--------------------------------------------------------------------------------
        ! Compute first moment
        !--------------------------------------------------------------------------------


        associate( &
            M  => first_moment, &
            f1 => integrant(:, :, 1), &
            f2 => integrant(:, :, 2), &
            f3 => integrant(:, :, 3) &
            )

            M%x = this%Compute_surface_integral( f1 )

            M%y = this%Compute_surface_integral( f2 )

            M%z = this%Compute_surface_integral( f3 )

        end associate

        !--------------------------------------------------------------------------------
        ! Deallocate array
        !--------------------------------------------------------------------------------

        deallocate ( &
            integrant, &
            stat   = deallocate_status, &
            errmsg = error_message )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (sphere_t)'
            write( stderr, '(A)' ) 'Deallocation failed in COMPUTE_FIRST_MOMENT'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end subroutine Compute_first_moment
    !
    !*****************************************************************************************
    !
    ! Public SPHEREPACK 3.2 methods
    !
    !*****************************************************************************************
    !
    subroutine Get_colatitude_derivative( this, polar_component, azimuthal_component )
        !
        !< Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vtsgs.html
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)   :: this
        real (WP),        intent (out)      :: polar_component(:)     !! vt
        real (WP),        intent (out)      :: azimuthal_component(:) !! wt
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        !        ! TODO incorporate Vtsgsi into type(workspace)
        !        subroutine vtsgsi(nlat, nlon, wvts, lwvts, work, lwork, dwork, ldwork, ierror)

        !        call Vtsgs( &
        !            size( this%grid%latitudes ), size( this%grid%longitudes ),  this%ityp, 1, &
        !            polar_component, azimuthal_component, &
        !            this%NLAT, this%NLON, &
        !            this%workspace%real_polar_harmonic_coefficients, this%workspace%imaginary_polar_harmonic_coefficients, &
        !            this%workspace%real_azimuthal_harmonic_coefficients, this%workspace%imaginary_azimuthal_harmonic_coefficients, &
        !            size(this%workspace%real_polar_harmonic_coefficients, dim = 1), size(this%workspace%real_polar_harmonic_coefficients, dim = 2), &
        !            this%workspace%wvts, this%workspace%size(wvts), &
        !            this%workspace%work, size(this%workspace%work), ierror)

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Get_colatitude_derivative
    !
    !*****************************************************************************************
    !
    subroutine Get_gradient( this, scalar_function, &
        polar_gradient_component, azimuthal_gradient_component )
        !
        ! SPHEREPACK 3.2 documentation
        !
        !     subroutine gradgs(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab,
        !                    wvhsgs, lvhsgs, work, lwork, ierror)
        !
        !     given the scalar spherical harmonic coefficients a and b, precomputed
        !     by subroutine shags for a scalar field sf, subroutine gradgs computes
        !     an irrotational vector field (v, w) such that
        !
        !           gradient(sf) = (v, w).
        !
        !     v is the colatitudinal and w is the east longitudinal component
        !     of the gradient.  i.e.,
        !
        !            v(i, j) = d(sf(i, j))/dtheta
        !
        !     and
        !
        !            w(i, j) = 1/sint*d(sf(i, j))/dlambda
        !
        !     at the gaussian colatitude point theta(i) (see nlat as input
        !     parameter) and longitude lambda(j) = (j-1)*2*pi/nlon where
        !     sint = sin(theta(i)).
        !
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are computed
        !            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid point
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater than
        !            3.  the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !
        !     isym   this has the same value as the isym that was input to
        !            subroutine shags to compute the arrays a and b from the
        !            scalar field sf.  isym determines whether (v, w) are
        !            computed on the full or half sphere as follows:
        !
        !      = 0
        !
        !           sf is not symmetric about the equator. in this case
        !           the vector field (v, w) is computed on the entire sphere.
        !           i.e., in the arrays  v(i, j), w(i, j) for i=1, ..., nlat and
        !           j=1, ..., nlon.
        !
        !      = 1
        !
        !           sf is antisymmetric about the equator. in this case w is
        !           antisymmetric and v is symmetric about the equator. w
        !           and v are computed on the northern hemisphere only.  i.e.,
        !           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
        !           and j=1, ..., nlon.  if nlat is even they are computed for
        !           i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !      = 2
        !
        !           sf is symmetric about the equator. in this case w is
        !           symmetric and v is antisymmetric about the equator. w
        !           and v are computed on the northern hemisphere only.  i.e.,
        !           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
        !           and j=1, ..., nlon.  if nlat is even they are computed for
        !           i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !
        !     nt     nt is the number of scalar and vector fields.  some
        !            computational efficiency is obtained for multiple fields.
        !            the arrays a, b, v, and w can be three dimensional corresponding
        !            to an indexed multiple array sf.  in this case, multiple
        !            vector synthesis will be performed to compute each vector
        !            field.  the third index for a, b, v, and w is the synthesis
        !            index which assumes the values k = 1, ..., nt.  for a single
        !            synthesis set nt = 1.  the description of the remaining
        !            parameters is simplified by assuming that nt=1 or that a, b, v,
        !            and w are two dimensional arrays.
        !
        !     idvw   the first dimension of the arrays v, w as it appears in
        !            the program that calls gradgs. if isym = 0 then idvw
        !            must be at least nlat.  if isym = 1 or 2 and nlat is
        !            even then idvw must be at least nlat/2. if isym = 1 or 2
        !            and nlat is odd then idvw must be at least (nlat+1)/2.
        !
        !     jdvw   the second dimension of the arrays v, w as it appears in
        !            the program that calls gradgs. jdvw must be at least nlon.
        !
        !     a, b    two or three dimensional arrays (see input parameter nt)
        !            that contain scalar spherical harmonic coefficients
        !            of the scalar field array sf as computed by subroutine shags.
        !     ***    a, b must be computed by shags prior to calling gradgs.
        !
        !     mdab   the first dimension of the arrays a and b as it appears in
        !            the program that calls gradgs (and shags). mdab must be at
        !            least min0(nlat, (nlon+2)/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     ndab   the second dimension of the arrays a and b as it appears in
        !            the program that calls gradgs (and shags). ndab must be at
        !            least nlat.
        !
        !
        !     wvhsgs an array which must be initialized by subroutine vhsgsi.
        !            once initialized,
        !            wvhsgs can be used repeatedly by gradgs as long as nlon
        !            and nlat remain unchanged.  wvhsgs must not be altered
        !            between calls of gradgs.
        !
        !
        !     lvhsgs the dimension of the array wvhsgs as it appears in the
        !            program that calls grradgs.  define
        !
        !               l1 = min0(nlat, nlon/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lvhsgs must be at least
        !
        !                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat
        !
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls gradgs. define
        !
        !               l1 = min0(nlat, nlon/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2                  if nlat is even or
        !               l2 = (nlat+1)/2              if nlat is odd
        !
        !            if isym = 0, lwork must be greater than or equal to
        !
        !               nlat*((2*nt+1)*nlon+2*l1*nt+1).
        !
        !            if isym = 1 or 2, lwork must be greater than or equal to
        !
        !               (2*nt+1)*l2*nlon+nlat*(2*l1*nt+1).
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !
        !     v, w   two or three dimensional arrays (see input parameter nt) that
        !           contain an irrotational vector field such that the gradient of
        !           the scalar field sf is (v, w).  w(i, j) is the east longitude
        !           component and v(i, j) is the colatitudinal component of velocity
        !           at gaussian colatitude and longitude lambda(j) = (j-1)*2*pi/nlon
        !           the indices for v and w are defined at the input parameter
        !           isym.  the vorticity of (v, w) is zero.  note that any nonzero
        !           vector field on the sphere will be multiple valued at the poles
        !           [reference swarztrauber].
        !
        !
        !  ierror   = 0  no errors
        !           = 1  error in the specification of nlat
        !           = 2  error in the specification of nlon
        !           = 3  error in the specification of isym
        !           = 4  error in the specification of nt
        !           = 5  error in the specification of idvw
        !           = 6  error in the specification of jdvw
        !           = 7  error in the specification of mdab
        !           = 8  error in the specification of ndab
        !           = 9  error in the specification of lvhsgs
        !           = 10 error in the specification of lwork
        !
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (in)      :: scalar_function(:, :)
        real (WP),        intent (out)     :: polar_gradient_component(:, :)
        real (WP),        intent (out)     :: azimuthal_gradient_component(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: error_flag
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Compute the (real) spherical harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%Perform_scalar_analysis( scalar_function )

        !--------------------------------------------------------------------------------
        ! Compute gradient
        !--------------------------------------------------------------------------------

        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            isym   => this%SCALAR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            v      => polar_gradient_component, &
            w      => azimuthal_gradient_component, &
            idvw   => this%NLAT, &
            jdvw   => this%NLON, &
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients, &
            mdab   => this%NLAT, &
            ndab   => this%NLAT, &
            wvhsgs => this%workspace%wvhsgs, &
            lvhsgs => size( this%workspace%wvhsgs ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            call Gradgs( nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
                wvhsgs, lvhsgs, work, lwork, ierror)

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: GET_GRADIENT'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of SCALAR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_SYNTHESES'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'POLAR_GRADIENT_COMPONENT (THETA)'
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'AZIMUTHAL_GRADIENT_COMPONENT (PHI)'
                write( stderr, '(A)') 'size( THETA, dim = 1)'
                write( stderr, '(A)') 'size( PHI,   dim = 1)'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'POLAR_GRADIENT_COMPONENT (THETA)'
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'AZIMUTHAL_GRADIENT_COMPONENT (PHI)'
                write( stderr, '(A)') 'size( THETA, dim = 2)'
                write( stderr, '(A)') 'size( PHI,   dim = 2)'


            else if ( error_flag == 7 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A)'
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
                write( stderr, '(A)') 'size( A, dim = 1)'
                write( stderr, '(A)') 'size( B, dim = 1)'


            else if (error_flag == 8) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A)'
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
                write( stderr, '(A)') 'size( A, dim = 2)'
                write( stderr, '(A)') 'size( B, dim = 2)'


            else if (error_flag == 9) then

                write( stderr, '(A)') 'Invalid extent for WVHSGS'
                write( stderr, '(A)') 'size( WVHSGS )'

            else if (error_flag == 10) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Get_gradient
    !
    !*****************************************************************************************
    !
    subroutine Invert_gradient( this, &
        polar_gradient_component, azimuthal_gradient_component, &
        scalar_function )
        !
        !     Documentation: SPHEREPACK 3.2
        !
        !     subroutine igradgs(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb,
        !                       wshsgs, lshsgs, work, lwork, ierror)
        !
        !     let br, bi, cr, ci be the vector spherical harmonic coefficients
        !     precomputed by vhags for a vector field (v, w).  let (v', w') be
        !     the irrotational component of (v, w) (i.e., (v', w') is generated
        !     by assuming cr, ci are zero and synthesizing br, bi with vhsgs).
        !     then subroutine igradgs computes a scalar field sf such that
        !
        !            gradient(sf) = (v', w').
        !
        !     i.e.,
        !
        !            v'(i, j) = d(sf(i, j))/dtheta          (colatitudinal component of
        !                                                 the gradient)
        !     and
        !
        !            w'(i, j) = 1/sint*d(sf(i, j))/dlambda  (east longitudinal component
        !                                                 of the gradient)
        !
        !     at the gaussian colatitude theta(i) (see nlat as input parameter)
        !     and longitude lambda(j) = (j-1)*2*pi/nlon where sint = sin(theta(i)).
        !
        !     note:  for an irrotational vector field (v, w), subroutine igradgs
        !     computes a scalar field whose gradient is (v, w).  in ay case,
        !     subroutine igradgs "inverts" the gradient subroutine gradgs.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are computed
        !            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid point
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater than
        !            3.  the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !
        !     isym   a parameter which determines whether the scalar field sf is
        !            computed on the full or half sphere as follows:
        !
        !      = 0
        !
        !            the symmetries/antsymmetries described in isym=1, 2 below
        !            do not exist in (v, w) about the equator.  in this case sf
        !            is neither symmetric nor antisymmetric about the equator.
        !            sf is computed on the entire sphere.  i.e., in the array
        !            sf(i, j) for i=1, ..., nlat and  j=1, ..., nlon
        !
        !      = 1
        !
        !            w is antisymmetric and v is symmetric about the equator.
        !            in this case sf is antisymmetyric about the equator and
        !            is computed for the northern hemisphere only.  i.e.,
        !            if nlat is odd sf is computed in the array sf(i, j) for
        !            i=1, ..., (nlat+1)/2 and for j=1, ..., nlon.  if nlat is even
        !            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
        !            and j=1, ..., nlon.
        !
        !      = 2
        !
        !            w is symmetric and v is antisymmetric about the equator.
        !            in this case sf is symmetyric about the equator and
        !            is computed for the northern hemisphere only.  i.e.,
        !            if nlat is odd sf is computed in the array sf(i, j) for
        !            i=1, ..., (nlat+1)/2 and for j=1, ..., nlon.  if nlat is even
        !            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
        !            and j=1, ..., nlon.
        !
        !
        !     nt     nt is the number of scalar and vector fields.  some
        !            computational efficiency is obtained for multiple fields.
        !            the arrays br, bi, and sf can be three dimensional corresponding
        !            to an indexed multiple vector field (v, w).  in this case,
        !            multiple scalar synthesis will be performed to compute each
        !            scalar field.  the third index for br, bi, and sf is the synthesis
        !            index which assumes the values k = 1, ..., nt.  for a single
        !            synthesis set nt = 1.  the description of the remaining
        !            parameters is simplified by assuming that nt=1 or that br, bi,
        !            and sf are two dimensional arrays.
        !
        !     isf    the first dimension of the array sf as it appears in
        !            the program that calls igradgs. if isym = 0 then isf
        !            must be at least nlat.  if isym = 1 or 2 and nlat is
        !            even then isf must be at least nlat/2. if isym = 1 or 2
        !            and nlat is odd then isf must be at least (nlat+1)/2.
        !
        !     jsf    the second dimension of the array sf as it appears in
        !            the program that calls igradgs. jsf must be at least nlon.
        !
        !     br, bi  two or three dimensional arrays (see input parameter nt)
        !            that contain vector spherical harmonic coefficients
        !            of the vector field (v, w) as computed by subroutine vhags.
        !     ***    br, bi must be computed by vhags prior to calling igradgs.
        !
        !     mdb    the first dimension of the arrays br and bi as it appears in
        !            the program that calls igradgs (and vhags). mdb must be at
        !            least min0(nlat, nlon/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     ndb    the second dimension of the arrays br and bi as it appears in
        !            the program that calls igradgs (and vhags). ndb must be at
        !            least nlat.
        !
        !
        !  wshsgs    an array which must be initialized by subroutine igradgsi
        !            (or equivalently by subroutine shsesi).  once initialized,
        !            wshsgs can be used repeatedly by igradgs as long as nlon
        !            and nlat remain unchanged.  wshsgs must not be altered
        !            between calls of igradgs.
        !
        !
        !  lshsgs    the dimension of the array wshsgs as it appears in the
        !            program that calls igradgs. define
        !
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd.
        !
        !
        !            then lshsgs must be greater than or equal to
        !
        !               nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
        !
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls igradgs. define
        !
        !               l2 = nlat/2                    if nlat is even or
        !               l2 = (nlat+1)/2                if nlat is odd
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            if isym = 0 lwork must be greater than or equal to
        !
        !               nlat*((nt+1)*nlon+2*nt*l1+1)
        !
        !            if isym > 0 lwork must be greater than or equal to
        !
        !               (nt+1)*l2*nlon+nlat*(2*nt*l1+1)
        !
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !
        !     sf    a two or three dimensional array (see input parameter nt) that
        !           contain a scalar field whose gradient is the irrotational
        !           component of the vector field (v, w).  the vector spherical
        !           harmonic coefficients br, bi were precomputed by subroutine
        !           vhags.  sf(i, j) is given at the gaussian colatitude theta(i)
        !           and longitude lambda(j) = (j-1)*2*pi/nlon.  the index ranges
        !           are defined at input parameter isym.
        !
        !
        !  ierror   = 0  no errors
        !           = 1  error in the specification of nlat
        !           = 2  error in the specification of nlon
        !           = 3  error in the specification of isym
        !           = 4  error in the specification of nt
        !           = 5  error in the specification of isf
        !           = 6  error in the specification of jsf
        !           = 7  error in the specification of mdb
        !           = 8  error in the specification of ndb
        !           = 9  error in the specification of lshsgs
        !           = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        real (WP),        intent (in)     :: polar_gradient_component(:, :)
        real (WP),        intent (in)     :: azimuthal_gradient_component(:, :)
        real (WP),        intent (out)    :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                      :: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Perform vector analysis
        !--------------------------------------------------------------------------------

         !!!!TODO

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

    !        associate( &
    !            nlat   => this%NLAT, &
    !            nlon   => this%NLON, &
    !            isym   => this%SCALAR_SYMMETRIES, &
    !            nt     => this%NUMBER_OF_SYNTHESES, &
    !            sf     => scalar_function, &
    !            isf    => size( scalar_function, dim = 1 ), &
    !            jsf    => size( scalar_function, dim = 2 ), &
    !            br     => this%real_polar_harmonic_coefficients, &
    !            bi     => this%imaginary_polar_harmonic_coefficients, &
    !            mdb    => size( this%real_polar_harmonic_coefficients, dim = 1 ), &
    !            ndb    => size( this%real_polar_harmonic_coefficients, dim = 2 ), &
    !            wshsgs => this%workspace%wshsgs, &
    !            lshsgs => size( this%workspace%wshsgs ), &
    !            work   => this%workspace%work, &
    !            lwork  => size( this%workspace%work ), &
    !            ierror => error_flag &
    !            )
    !
    !            call Igradgs( nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, &
    !                wshsgs, lshsgs, work, lwork, ierror )
    !
    !        end associate

        !--------------------------------------------------------------------------------
        ! Address error flag
        !--------------------------------------------------------------------------------

         !!!!TODO

    end subroutine Invert_gradient
    !
    !*****************************************************************************************
    !
    subroutine Get_divergence( this, vector_field, divergence )
        !
        ! Reference:
        ! http://www2.cisl.ucar.edu/spherepack/documentation#divgs.html
        !
        ! Documentation: SPHEREPACK 3.2
        !
        !
        !     subroutine divgs(nlat, nlon, isym, nt, divg, idiv, jdiv, br, bi, mdb, ndb,
        !                    wshsgs, lshsgs, work, lwork, ierror)
        !
        !     given the vector spherical harmonic coefficients br and bi, precomputed
        !     by subroutine vhags for a vector field (v, w), subroutine divgs
        !     computes the divergence of the vector field in the scalar array divg.
        !     divg(i, j) is the divergence at the gaussian colatitude point theta(i)
        !     (see nlat as input parameter) and east longitude
        !
        !            lambda(j) = (j-1)*2*pi/nlon
        !
        !     on the sphere.  i.e.
        !
        !            dv(i, j) = 1/sint*[ d(sint*v(i, j))/dtheta + d(w(i, j))/dlambda ]
        !
        !     where sint = sin(theta(i)).  w is the east longitudinal and v
        !     is the colatitudinal component of the vector field from which
        !     br, bi were precomputed
        !
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are computed
        !            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid point
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than zero. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !
        !     isym   a parameter which determines whether the divergence is
        !            computed on the full or half sphere as follows:
        !
        !      = 0
        !
        !            the symmetries/antsymmetries described in isym=1, 2 below
        !            do not exist in (v, w) about the equator.  in this case the
        !            divergence is neither symmetric nor antisymmetric about
        !            the equator.  the divergence is computed on the entire
        !            sphere.  i.e., in the array divg(i, j) for i=1, ..., nlat and
        !            j=1, ..., nlon.
        !
        !      = 1
        !
        !            w is antisymmetric and v is symmetric about the equator.
        !            in this case the divergence is antisymmetyric about
        !            the equator and is computed for the northern hemisphere
        !            only.  i.e., if nlat is odd the divergence is computed
        !            in the array divg(i, j) for i=1, ..., (nlat+1)/2 and for
        !            j=1, ..., nlon.  if nlat is even the divergence is computed
        !            in the array divg(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !      = 2
        !            w is symmetric and v is antisymmetric about the equator
        !            in this case the divergence is symmetyric about the
        !            equator and is computed for the northern hemisphere
        !            only.  i.e., if nlat is odd the divergence is computed
        !            in the array divg(i, j) for i=1, ..., (nlat+1)/2 and for
        !            j=1, ..., nlon.  if nlat is even the divergence is computed
        !            in the array divg(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !
        !     nt     nt is the number of scalar and vector fields.  some
        !            computational efficiency is obtained for multiple fields.
        !            in the program that calls divgs, the arrays br, bi, and divg
        !            can be three dimensional corresponding to an indexed multiple
        !            vector field.  in this case multiple scalar synthesis will
        !            be performed to compute the divergence for each field.  the
        !            third index is the synthesis index which assumes the values
        !            k=1, ..., nt.  for a single synthesis set nt = 1.  the
        !            description of the remaining parameters is simplified by
        !            assuming that nt=1 or that all the arrays are two dimensional.
        !
        !     idiv   the first dimension of the array divg as it appears in
        !            the program that calls divgs. if isym = 0 then idiv
        !            must be at least nlat.  if isym = 1 or 2 and nlat is
        !            even then idiv must be at least nlat/2. if isym = 1 or 2
        !            and nlat is odd then idiv must be at least (nlat+1)/2.
        !
        !     jdiv   the second dimension of the array divg as it appears in
        !            the program that calls divgs. jdiv must be at least nlon.
        !
        !     br, bi  two or three dimensional arrays (see input parameter nt)
        !            that contain vector spherical harmonic coefficients
        !            of the vector field (v, w) as computed by subroutine vhags.
        !     ***    br and bi must be computed by vhags prior to calling
        !            divgs.
        !
        !     mdb    the first dimension of the arrays br and bi as it
        !            appears in the program that calls divgs. mdb must be at
        !            least min0(nlat, nlon/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     ndb    the second dimension of the arrays br and bi as it
        !            appears in the program that calls divgs. ndb must be at
        !            least nlat.
        !
        !
        !     wshsgs an array which must be intialized by subroutine shsgsi.
        !            once initialized,
        !            wshsgs can be used repeatedly by divgs as long as nlon
        !            and nlat remain unchanged.  wshsgs must not be altered
        !            between calls of divgs.
        !
        !
        !     lshsgs the dimension of the array wshsgs as it appears in the
        !            program that calls divgs. define
        !
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lshsgs must be at least
        !
        !               nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
        !
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls divgs. define
        !
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2                    if nlat is even or
        !               l2 = (nlat+1)/2                if nlat is odd
        !
        !            if isym = 0 then lwork must be at least
        !
        !               nlat*((nt+1)*nlon+2*nt*l1+1)
        !
        !            if isym > 0 then lwork must be at least
        !
        !               (nt+1)*l2*nlon+nlat*(2*nt*l1+1)
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !
        !    divg   a two or three dimensional array (see input parameter nt)
        !           that contains the divergence of the vector field (v, w)
        !           whose coefficients br, bi where computed by subroutine
        !           vhags.  divg(i, j) is the divergence at the gaussian colatitude
        !           point theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon.
        !           the index ranges are defined above at the input parameter
        !           isym.
        !
        !
        !    ierror = 0  no errors
        !           = 1  error in the specification of nlat
        !           = 2  error in the specification of nlon
        !           = 3  error in the specification of isym
        !           = 4  error in the specification of nt
        !           = 5  error in the specification of idiv
        !           = 6  error in the specification of jdiv
        !           = 7  error in the specification of mdb
        !           = 8  error in the specification of ndb
        !           = 9  error in the specification of lshsgs
        !           = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        real (WP),        intent (in)     :: vector_field (:, :, :)
        real (WP),        intent (out)    :: divergence (:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Calculate the (real) vector harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%Perform_vector_analysis( vector_field )

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            isym   => this%SCALAR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            divg   => divergence, &
            idiv   => size( divergence, dim = 1), &
            jdiv   => size( divergence, dim = 2), &
            br     => this%workspace%real_polar_harmonic_coefficients, &
            bi     => this%workspace%imaginary_polar_harmonic_coefficients, &
            mdb    => size( this%workspace%real_polar_harmonic_coefficients, dim = 1 ), &
            ndb    => size( this%workspace%real_polar_harmonic_coefficients, dim = 2 ), &
            wshsgs => this%workspace%wshsgs, &
            lshsgs => size( this%workspace%wshsgs ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            call Divgs( nlat, nlon, isym, nt, divg, idiv, jdiv, br, bi, mdb, ndb, &
                wshsgs, lshsgs, work, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: GET_DIVERGENCE'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of SCALAR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_SYNTHESES'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for DIVERGENCE'
                write( stderr, '(A)') 'size( DIVERGENCE, dim = 1)'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Invalid extent for DIVERGENCE'
                write( stderr, '(A)') 'size( DIVERGENCE, dim = 2)'

            else if ( error_flag == 7 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_POLAR_HARMONIC_COEFFICIENTS (BR)'
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_POLAR_HARMONIC_COEFFICIENTS (BI)'
                write( stderr, '(A)') 'size( BR, dim = 1)'
                write( stderr, '(A)') 'size( BI, dim = 1)'


            else if (error_flag == 8) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_POLAR_HARMONIC_COEFFICIENTS (BR)'
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_POLAR_HARMONIC_COEFFICIENTS (BI)'
                write( stderr, '(A)') 'size( BR, dim = 2)'
                write( stderr, '(A)') 'size( BI, dim = 2)'


            else if (error_flag == 9) then

                write( stderr, '(A)') 'Invalid extent for WSHSGS'
                write( stderr, '(A)') 'size( WSHSGS )'

            else if (error_flag == 10) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if


    end subroutine Get_divergence
    !
    !*****************************************************************************************
    !
    subroutine Invert_divergence( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Invert_divergence
    !
    !*****************************************************************************************
    !
    subroutine Get_vorticity( this, vector_function, vorticity )
        !
        !--------------------------------------------------------------------------------
        !
        !    Documentation: SPHEREPACK 3.2
        !
        !    subroutine vrtgs(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc,
        !                     wshsgs, lshsgs, work, lwork, ierror)
        !
        !     given the vector spherical harmonic coefficients cr and ci, precomputed
        !     by subroutine vhags for a vector field (v, w), subroutine vrtgs
        !     computes the vorticity of the vector field in the scalar array
        !     vort.  vort(i, j) is the vorticity at the gaussian colatitude
        !     theta(i) (see nlat as input parameter) and longitude
        !     lambda(j) = (j-1)*2*pi/nlon on the sphere.  i.e.,
        !
        !            vort(i, j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint
        !
        !     where sint = sin(theta(i)).  w is the east longitudinal and v
        !     is the colatitudinal component of the vector field from which
        !     cr, ci were precomputed.  required associated legendre polynomials
        !     are stored rather than recomputed as they are in subroutine vrtgc.
        !
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are computed
        !            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid point
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than 3. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !
        !     isym   a parameter which determines whether the vorticity is
        !            computed on the full or half sphere as follows:
        !
        !      = 0
        !            the symmetries/antsymmetries described in isym=1, 2 below
        !            do not exist in (v, w) about the equator.  in this case the
        !            vorticity is neither symmetric nor antisymmetric about
        !            the equator.  the vorticity is computed on the entire
        !            sphere.  i.e., in the array vort(i, j) for i=1, ..., nlat and
        !            j=1, ..., nlon.
        !
        !      = 1
        !            w is antisymmetric and v is symmetric about the equator.
        !            in this case the vorticity is symmetyric about the
        !            equator and is computed for the northern hemisphere
        !            only.  i.e., if nlat is odd the vorticity is computed
        !            in the array vort(i, j) for i=1, ..., (nlat+1)/2 and for
        !            j=1, ..., nlon.  if nlat is even the vorticity is computed
        !            in the array vort(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !      = 2
        !            w is symmetric and v is antisymmetric about the equator
        !            in this case the vorticity is antisymmetric about the
        !            equator and is computed for the northern hemisphere
        !            only.  i.e., if nlat is odd the vorticity is computed
        !            in the array vort(i, j) for i=1, ..., (nlat+1)/2 and for
        !            j=1, ..., nlon.  if nlat is even the vorticity is computed
        !            in the array vort(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !
        !      nt    nt is the number of scalar and vector fields.  some
        !            computational efficiency is obtained for multiple fields.
        !            in the program that calls vrtgs, the arrays cr, ci, and vort
        !            can be three dimensional corresponding to an indexed multiple
        !            vector field.  in this case multiple scalar synthesis will
        !            be performed to compute the vorticity for each field.  the
        !            third index is the synthesis index which assumes the values
        !            k=1, ..., nt.  for a single synthesis set nt = 1.  the
        !            description of the remaining parameters is simplified by
        !            assuming that nt=1 or that all the arrays are two dimensional.
        !
        !     ivrt   the first dimension of the array vort as it appears in
        !            the program that calls vrtgs. if isym = 0 then ivrt
        !            must be at least nlat.  if isym = 1 or 2 and nlat is
        !            even then ivrt must be at least nlat/2. if isym = 1 or 2
        !            and nlat is odd then ivrt must be at least (nlat+1)/2.
        !
        !     jvrt   the second dimension of the array vort as it appears in
        !            the program that calls vrtgs. jvrt must be at least nlon.
        !
        !    cr, ci   two or three dimensional arrays (see input parameter nt)
        !            that contain vector spherical harmonic coefficients
        !            of the vector field (v, w) as computed by subroutine vhags.
        !     ***    cr and ci must be computed by vhags prior to calling
        !            vrtgs.
        !
        !      mdc   the first dimension of the arrays cr and ci as it
        !            appears in the program that calls vrtgs. mdc must be at
        !            least min0(nlat, nlon/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !      ndc   the second dimension of the arrays cr and ci as it
        !            appears in the program that calls vrtgs. ndc must be at
        !            least nlat.
        !
        !   wshsgs   an array which must be initialized by subroutine shsgsi.
        !            once initialized,
        !            wshsgs can be used repeatedly by vrtgs as long as nlon
        !            and nlat remain unchanged.  wshsgs must not be altered
        !            between calls of vrtgs
        !
        !   lshsgs   the dimension of the array wshsgs   as it appears in the
        !            program that calls vrtgs. define
        !
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lshsgs must be at least
        !
        !            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
        !
        !     work   a work array that does not have to be saved.
        !
        !    lwork   the dimension of the array work as it appears in the
        !            program that calls vrtgs. define
        !
        !               l1 = min0(nlat, nlon/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd.
        !
        !            if isym = 0 then lwork must be at least
        !
        !               nlat*((nt+1)*nlon+2*nt*l1+1)
        !
        !            if isym > 0 then lwork must be at least
        !
        !               (nt+1)*l2*nlon+nlat*(2*nt*l1+1)
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !
        !     vort   a two or three dimensional array (see input parameter nt)
        !            that contains the vorticity of the vector field (v, w)
        !            whose coefficients cr, ci where computed by subroutine vhags.
        !            vort(i, j) is the vorticity at the gaussian colatitude point
        !            theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon.
        !            the index ranges are defined above at the input parameter
        !            isym.
        !
        !
        !   ierror   an error parameter which indicates fatal errors with input
        !            parameters when returned positive.
        !          = 0  no errors
        !          = 1  error in the specification of nlat
        !          = 2  error in the specification of nlon
        !          = 3  error in the specification of isym
        !          = 4  error in the specification of nt
        !          = 5  error in the specification of ivrt
        !          = 6  error in the specification of jvrt
        !          = 7  error in the specification of mdc
        !          = 8  error in the specification of ndc
        !          = 9  error in the specification of lshsgs
        !          = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (in)      :: vector_function(:, :, :)
        real (WP),        intent (out)     :: vorticity(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                        :: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! calculate the (real) vector harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%Perform_vector_analysis( vector_function )

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2 routine
        !--------------------------------------------------------------------------------

        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            isym   => this%SCALAR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            vort   => vorticity, &
            ivrt   => size( vorticity, dim = 1 ), &
            jvrt   => size( vorticity, dim = 2 ), &
            cr     => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci     => this%workspace%imaginary_azimuthal_harmonic_coefficients, &
            mdc    => size( this%workspace%real_azimuthal_harmonic_coefficients, dim = 1 ), &
            ndc    => size( this%workspace%real_azimuthal_harmonic_coefficients, dim = 2 ), &
            wshsgs => this%workspace%wshsgs, &
            lshsgs => size( this%workspace%wshsgs ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            ! calculate the surface vorticity
            call Vrtgs( nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
                wshsgs, lshsgs, work, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'ERROR: GET_VORTICITY'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'


            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of SCALAR_SYMMETRIES'


            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_SYNTHESES'



            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for SOURCE_TERM'
                write( stderr, '(A)') 'size( SOURCE_TERM, dim = 1)'


            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Invalid extent for SOURCE_TERM'
                write( stderr, '(A)') 'size( SOURCE_TERM, dim = 2)'


            else if ( error_flag == 7 ) then

                write( stderr, '(A)') 'Invalid extent for CR'
                write( stderr, '(A)') 'size( CR, dim = 1)'


            else if (error_flag == 8) then

                write( stderr, '(A)') 'Invalid extent for CR'
                write( stderr, '(A)') 'size( CR, dim = 2)'


            else if (error_flag == 9) then

                write( stderr, '(A)') 'Invalid extent for WSHSGS'
                write( stderr, '(A)') 'size( WSHSGS )'

            else if (error_flag == 10) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Get_vorticity
    !
    !*****************************************************************************************
    !
    subroutine Invert_vorticity( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Invert_vorticity
    !
    !*****************************************************************************************
    !
    subroutine Invert_divergence_and_vorticity( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Invert_divergence_and_vorticity
    !
    !*****************************************************************************************
    !
    subroutine Get_scalar_laplacian( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Get_scalar_laplacian
    !
    !*****************************************************************************************
    !
    subroutine Invert_helmholtz( this, helmholtz_constant, source_term, solution )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)         :: this
        real (WP), intent (in)                    :: helmholtz_constant
        real (WP), dimension (:, :), intent (in)   :: source_term
        real (WP), dimension (:, :), intent (out)  :: solution
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (WP)     :: perturbation
        integer (IP)  :: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        call this%Perform_scalar_analysis( source_term )

        ! invert the helmholtz (or poisson) equation
        !
        ! see: https://www2.cisl.ucar.edu/spherepack/documentation#islapgs.html
        call Islapgs( &
            this%NLAT, this%NLON, &
            this%SCALAR_SYMMETRIES, 1, helmholtz_constant, &
            solution, this%NLAT, this%NLON, &
            this%workspace%real_harmonic_coefficients, this%workspace%imaginary_harmonic_coefficients, &
            this%NLAT, this%NLAT, &
            this%workspace%wshsgs, size( this%workspace%wshsgs ), &
            this%workspace%work, size( this%workspace%work ), &
            perturbation, ierror)

        ! check for errors
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Islapgs'
            return
        end if

    end subroutine Invert_helmholtz
    !
    !*****************************************************************************************
    !
    subroutine Get_vector_laplacian( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Get_vector_laplacian
    !
    !*****************************************************************************************
    !
    subroutine Invert_vector_laplacian( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Invert_vector_laplacian
    !
    !*****************************************************************************************
    !
    subroutine Get_stream_function_and_velocity_potential( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Get_stream_function_and_velocity_potential
    !
    !*****************************************************************************************
    !
    subroutine Invert_stream_function_and_velocity_potential( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Invert_stream_function_and_velocity_potential
    !
    !*****************************************************************************************
    !
    subroutine Perform_grid_transfers( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Perform_grid_transfers
    !
    !*****************************************************************************************
    !
    subroutine Perform_geo_math_coordinate_transfers( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Perform_geo_math_coordinate_transfers
    !
    !*****************************************************************************************
    !
    subroutine Perform_scalar_analysis( this, scalar_function )
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        !
        ! SPHEREPACK 3.2 documentation
        !
        !     subroutine shags(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
        !                      wshags, lshags, work, lwork, ierror)
        !
        !     subroutine shags performs the spherical harmonic analysis
        !     on the array g and stores the result in the arrays a and b.
        !     the analysis is performed on a gaussian grid in colatitude
        !     and an equally spaced grid in longitude.  the associated
        !     legendre functions are stored rather than recomputed as they
        !     are in subroutine shagc.  the analysis is described below
        !     at output parameters a, b.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are compu
        !            in radians in theta(1), ..., theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid poi
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than or equal to 4. the efficiency of the computation is
        !            improved when nlon is a product of small prime numbers.
        !
        !     isym   = 0  no symmetries exist about the equator. the analysis
        !                 is performed on the entire sphere.  i.e. on the
        !                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
        !                 (see description of g below)
        !
        !            = 1  g is antisymmetric about the equator. the analysis
        !                 is performed on the northern hemisphere only.  i.e.
        !                 if nlat is odd the analysis is performed on the
        !                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
        !                 if nlat is even the analysis is performed on the
        !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !
        !            = 2  g is symmetric about the equator. the analysis is
        !                 performed on the northern hemisphere only.  i.e.
        !                 if nlat is odd the analysis is performed on the
        !                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
        !                 if nlat is even the analysis is performed on the
        !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !     nt     the number of analyses.  in the program that calls shags,
        !            the arrays g, a and b can be three dimensional in which
        !            case multiple analyses will be performed.  the third
        !            index is the analysis index which assumes the values
        !            k=1, ..., nt.  for a single analysis set nt=1. the
        !            discription of the remaining parameters is simplified
        !            by assuming that nt=1 or that the arrays g, a and b
        !            have only two dimensions.
        !
        !     g      a two or three dimensional array (see input parameter
        !            nt) that contains the discrete function to be analyzed.
        !            g(i, j) contains the value of the function at the gaussian
        !            point theta(i) and longitude point phi(j) = (j-1)*2*pi/nlon
        !            the index ranges are defined above at the input parameter
        !            isym.
        !
        !     idg    the first dimension of the array g as it appears in the
        !            program that calls shags. if isym equals zero then idg
        !            must be at least nlat.  if isym is nonzero then idg must
        !            be at least nlat/2 if nlat is even or at least (nlat+1)/2
        !            if nlat is odd.
        !
        !     jdg    the second dimension of the array g as it appears in the
        !            program that calls shags. jdg must be at least nlon.
        !
        !     mdab   the first dimension of the arrays a and b as it appears
        !            in the program that calls shags. mdab must be at least
        !            min0((nlon+2)/2, nlat) if nlon is even or at least
        !            min0((nlon+1)/2, nlat) if nlon is odd.
        !
        !     ndab   the second dimension of the arrays a and b as it appears
        !            in the program that calls shags. ndab must be at least nlat
        !
        !     wshags an array which must be initialized by subroutine shagsi.
        !            once initialized, wshags can be used repeatedly by shags
        !            as long as nlat and nlon remain unchanged.  wshags must
        !            not be altered between calls of shags.
        !
        !     lshags the dimension of the array wshags as it appears in the
        !            program that calls shags. define
        !
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lshags must be at least
        !
        !            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
        !
        !     work   a real work space which need not be saved
        !
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls shags. define
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !
        !            if isym is zero then lwork must be at least
        !
        !                  nlat*nlon*(nt+1)
        !
        !            if isym is nonzero then lwork must be at least
        !
        !                  l2*nlon*(nt+1)
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     a, b    both a, b are two or three dimensional arrays (see input
        !            parameter nt) that contain the spherical harmonic
        !            coefficients in the representation of g(i, j) given in the
        !            discription of subroutine shags. for isym=0, a(m, n) and
        !            b(m, n) are given by the equations listed below. symmetric
        !            versions are used when isym is greater than zero.
        !
        !     definitions
        !
        !     1. the normalized associated legendre functions
        !
        !     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
        !                       *sin(theta)**m/(2**n*factorial(n)) times the
        !                       (n+m)th derivative of (x**2-1)**n with respect
        !                       to x=cos(theta).
        !
        !     2. the fourier transform of g(i, j).
        !
        !     c(m, i)          = 2/nlon times the sum from j=1 to j=nlon of
        !                       g(i, j)*cos((m-1)*(j-1)*2*pi/nlon)
        !                       (the first and last terms in this sum
        !                       are divided by 2)
        !
        !     s(m, i)          = 2/nlon times the sum from j=2 to j=nlon of
        !                       g(i, j)*sin((m-1)*(j-1)*2*pi/nlon)
        !
        !
        !     3. the gaussian points and weights on the sphere
        !        (computed by subroutine gaqd).
        !
        !        theta(1), ..., theta(nlat) (gaussian pts in radians)
        !        wts(1), ..., wts(nlat) (corresponding gaussian weights)
        !
        !
        !     4. the maximum (plus one) longitudinal wave number
        !
        !            mmax = min0(nlat, (nlon+2)/2) if nlon is even or
        !            mmax = min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !
        !     then for m=0, ..., mmax-1 and n=m, ..., nlat-1 the arrays a, b
        !     are given by
        !
        !     a(m+1, n+1)     =  the sum from i=1 to i=nlat of
        !                       c(m+1, i)*wts(i)*pbar(m, n, theta(i))
        !
        !     b(m+1, n+1)      = the sum from i=1 to nlat of
        !                       s(m+1, i)*wts(i)*pbar(m, n, theta(i))
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of isym
        !            = 4  error in the specification of nt
        !            = 5  error in the specification of idg
        !            = 6  error in the specification of jdg
        !            = 7  error in the specification of mdab
        !            = 8  error in the specification of ndab
        !            = 9  error in the specification of lshags
        !            = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)         :: this
        real (WP), dimension (:, :), intent (in)  :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: error_flag
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! perform the (real) spherical harmonic analysis
        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            isym   => this%SCALAR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            g      => scalar_function, &
            idg    => this%NLAT, &
            jdg    => this%NLON, &
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients, &
            mdab   => this%NLAT, &
            ndab   => this%NLAT, &
            wshags => this%workspace%wshags, &
            lshags => size( this%workspace%wshags ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            call Shags( nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                wshags, lshags, work, lwork, ierror )

        end associate

        ! check the error status
        if ( error_flag /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', error_flag, ' in Shags'
            return
        end if

    end subroutine Perform_scalar_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_scalar_synthesis( this, scalar_function )
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shsgs.html
        !
        ! SPHEREPACK 3.2 documentation
        !     subroutine shsgs(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
        !                        wshsgs, lshsgs, work, lwork, ierror)
        !
        !     subroutine shsgs performs the spherical harmonic synthesis
        !     on the arrays a and b and stores the result in the array g.
        !     the synthesis is performed on an equally spaced longitude grid
        !     and a gaussian colatitude grid.  the associated legendre functions
        !     are stored rather than recomputed as they are in subroutine
        !     shsgc.  the synthesis is described below at output parameter
        !     g.
        !
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are compu
        !            in radians in theta(1), ..., theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid poi
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than or equal to 4. the efficiency of the computation is
        !            improved when nlon is a product of small prime numbers.
        !
        !     isym   = 0  no symmetries exist about the equator. the synthesis
        !                 is performed on the entire sphere.  i.e. on the
        !                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
        !                 (see description of g below)
        !
        !            = 1  g is antisymmetric about the equator. the synthesis
        !                 is performed on the northern hemisphere only.  i.e.
        !                 if nlat is odd the synthesis is performed on the
        !                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
        !                 if nlat is even the synthesis is performed on the
        !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !
        !            = 2  g is symmetric about the equator. the synthesis is
        !                 performed on the northern hemisphere only.  i.e.
        !                 if nlat is odd the synthesis is performed on the
        !                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
        !                 if nlat is even the synthesis is performed on the
        !                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !     nt     the number of syntheses.  in the program that calls shsgs,
        !            the arrays g, a and b can be three dimensional in which
        !            case multiple synthesis will be performed.  the third
        !            index is the synthesis index which assumes the values
        !            k=1, ..., nt.  for a single synthesis set nt=1. the
        !            discription of the remaining parameters is simplified
        !            by assuming that nt=1 or that the arrays g, a and b
        !            have only two dimensions.
        !
        !     idg    the first dimension of the array g as it appears in the
        !            program that calls shagc. if isym equals zero then idg
        !            must be at least nlat.  if isym is nonzero then idg must
        !            be at least nlat/2 if nlat is even or at least (nlat+1)/2
        !            if nlat is odd.
        !
        !     jdg    the second dimension of the array g as it appears in the
        !            program that calls shagc. jdg must be at least nlon.
        !
        !     a, b    two or three dimensional arrays (see the input parameter
        !            nt) that contain the coefficients in the spherical harmonic
        !            expansion of g(i, j) given below at the definition of the
        !            output parameter g.  a(m, n) and b(m, n) are defined for
        !            indices m=1, ..., mmax and n=m, ..., nlat where mmax is the
        !            maximum (plus one) longitudinal wave number given by
        !            mmax = min0(nlat, (nlon+2)/2) if nlon is even or
        !            mmax = min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     mdab   the first dimension of the arrays a and b as it appears
        !            in the program that calls shsgs. mdab must be at least
        !            min0((nlon+2)/2, nlat) if nlon is even or at least
        !            min0((nlon+1)/2, nlat) if nlon is odd.
        !
        !     ndab   the second dimension of the arrays a and b as it appears
        !            in the program that calls shsgs. ndab must be at least nlat
        !
        !     wshsgs an array which must be initialized by subroutine shsgsi.
        !            once initialized, wshsgs can be used repeatedly by shsgs
        !            as long as nlat and nlon remain unchanged.  wshsgs must
        !            not be altered between calls of shsgs.
        !
        !     lshsgs the dimension of the array wshsgs as it appears in the
        !            program that calls shsgs. define
        !
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lshsgs must be at least
        !
        !            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
        !
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls shsgs. define
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !
        !            if isym is zero then lwork must be at least
        !
        !                  nlat*nlon*(nt+1)
        !
        !            if isym is nonzero then lwork must be at least
        !
        !                  l2*nlon*(nt+1)
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     g      a two or three dimensional array (see input parameter nt)
        !            that contains the discrete function which is synthesized.
        !            g(i, j) contains the value of the function at the gaussian
        !            colatitude point theta(i) and longitude point
        !            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined
        !            above at the input parameter isym.  for isym=0, g(i, j)
        !            is given by the the equations listed below.  symmetric
        !            versions are used when isym is greater than zero.
        !
        !     the normalized associated legendre functions are given by
        !
        !     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
        !                       *sin(theta)**m/(2**n*factorial(n)) times the
        !                       (n+m)th derivative of (x**2-1)**n with respect
        !                       to x=cos(theta)
        !
        !     define the maximum (plus one) longitudinal wave number
        !     as   mmax = min0(nlat, (nlon+2)/2) if nlon is even or
        !          mmax = min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     then g(i, j) = the sum from n=0 to n=nlat-1 of
        !
        !                   .5*pbar(0, n, theta(i))*a(1, n+1)
        !
        !              plus the sum from m=1 to m=mmax-1 of
        !
        !                   the sum from n=m to n=nlat-1 of
        !
        !              pbar(m, n, theta(i))*(a(m+1, n+1)*cos(m*phi(j))
        !                                    -b(m+1, n+1)*sin(m*phi(j)))
        !
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of isym
        !            = 4  error in the specification of nt
        !            = 5  error in the specification of idg
        !            = 6  error in the specification of jdg
        !            = 7  error in the specification of mdab
        !            = 8  error in the specification of ndab
        !            = 9  error in the specification of lshsgs
        !            = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (out)     :: scalar_function(:,:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: error_flag
        !--------------------------------------------------------------------------------

        ! perform (real) spherical harmonic synthesis
        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            isym   => this%SCALAR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            g      => scalar_function, &
            idg    => size( scalar_function, dim = 1), &
            jdg    => size( scalar_function, dim = 2), &
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients, &
            mdab   => size( this%workspace%real_harmonic_coefficients, dim = 1), &
            ndab   => size( this%workspace%real_harmonic_coefficients, dim = 2), &
            wshsgs => this%workspace%wshsgs, &
            lshsgs => size( this%workspace%wshsgs ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            call Shsgs( nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                wshsgs, lshsgs, work, lwork, ierror )

        end associate

        ! check the error status
        if ( error_flag /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', error_flag, ' in Shsgs'
            return
        end if

    end subroutine Perform_scalar_synthesis
    !
    !*****************************************************************************************
    ! TODO
    subroutine Perform_scalar_projection( this, scalar_function, scalar_projection )
        !
        !< Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shpg.html
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:, :), intent (in)  :: scalar_function
        real (WP), dimension (:, :), intent (out) :: scalar_projection
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! TODO: Include driver program into type(workspace)
        !        call Shpgi( &
        !            this%NLAT, this%NLON, this%SCALAR_SYMMETRIES, this%NTRUNC, &
        !            this%workspace%wshp, size( this%workspace%wshp ), &
        !            this%workspace%iwshp, size( this%workspace%iwshp ), &
        !            this%workspace%work, size( this%workspace%work ), ierror )
        !
        !        call Shpg( &
        !            this%NLAT, this%NLON, this%SCALAR_SYMMETRIES, this%NTRUNC, &
        !            scalar_function, scalar_projection, this%NLAT, &
        !            this%workspace%wshp, size( this%workspace%wshp ), &
        !            this%workspace%iwshp, size( this%workspace%iwshp ), &
        !            this%workspace%work, size( this%workspace%work ), ierror )

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Perform_scalar_projection
    !
    !*****************************************************************************************
    !
    subroutine Perform_vector_analysis( this, vector_function )
        !
        ! Reference:
        ! http://www2.cisl.ucar.edu/spherepack/documentation#vhags.html
        !
        ! Documentation: SPHEREPACK 3.2
        !
        !
        !     subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci,
        !                     mdab, ndab, wvhags, lvhags, work, lwork, ierror)
        !
        !     subroutine vhags performs the vector spherical harmonic analysis
        !     on the vector field (v, w) and stores the result in the arrays
        !     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
        !     (measured from the north pole) and east longitudinal components
        !     respectively, located at the gaussian colatitude point theta(i)
        !     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
        !     representation of (v, w) is given at output parameters v, w in
        !     subroutine vhses.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are computed
        !            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid point
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than zero. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !     ityp   = 0  no symmetries exist about the equator. the analysis
        !                 is performed on the entire sphere.  i.e. on the
        !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon.
        !
        !            = 1  no symmetries exist about the equator. the analysis
        !                 is performed on the entire sphere.  i.e. on the
        !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon. the curl of (v, w) is zero. that is,
        !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
        !                 the coefficients cr and ci are zero.
        !
        !            = 2  no symmetries exist about the equator. the analysis
        !                 is performed on the entire sphere.  i.e. on the
        !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e.,
        !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
        !                 the coefficients br and bi are zero.
        !
        !            = 3  v is symmetric and w is antisymmetric about the
        !                 equator. the analysis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the analysis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the analysis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !            = 4  v is symmetric and w is antisymmetric about the
        !                 equator. the analysis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the analysis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the analysis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the curl of (v, w) is zero. that is,
        !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
        !                 the coefficients cr and ci are zero.
        !
        !            = 5  v is symmetric and w is antisymmetric about the
        !                 equator. the analysis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the analysis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the analysis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the divergence of (v, w) is zero. i.e.,
        !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
        !                 the coefficients br and bi are zero.
        !
        !            = 6  v is antisymmetric and w is symmetric about the
        !                 equator. the analysis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the analysis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the analysis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !            = 7  v is antisymmetric and w is symmetric about the
        !                 equator. the analysis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the analysis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the analysis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the curl of (v, w) is zero. that is,
        !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
        !                 the coefficients cr and ci are zero.
        !
        !            = 8  v is antisymmetric and w is symmetric about the
        !                 equator. the analysis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the analysis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the analysis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the divergence of (v, w) is zero. i.e.,
        !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
        !                 the coefficients br and bi are zero.
        !
        !
        !     nt     the number of analyses.  in the program that calls vhags,
        !            the arrays v, w, br, bi, cr, and ci can be three dimensional
        !            in which case multiple analyses will be performed.
        !            the third index is the analysis index which assumes the
        !            values k=1, ..., nt.  for a single analysis set nt=1. the
        !            discription of the remaining parameters is simplified
        !            by assuming that nt=1 or that all the arrays are two
        !            dimensional.
        !
        !     v, w    two or three dimensional arrays (see input parameter nt)
        !            that contain the vector function to be analyzed.
        !            v is the colatitudnal component and w is the east
        !            longitudinal component. v(i, j), w(i, j) contain the
        !            components at the gaussian colatitude point theta(i)
        !            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
        !            are defined above at the input parameter ityp.
        !
        !     idvw   the first dimension of the arrays v, w as it appears in
        !            the program that calls vhags. if ityp .le. 2 then idvw
        !            must be at least nlat.  if ityp .gt. 2 and nlat is
        !            even then idvw must be at least nlat/2. if ityp .gt. 2
        !            and nlat is odd then idvw must be at least (nlat+1)/2.
        !
        !     jdvw   the second dimension of the arrays v, w as it appears in
        !            the program that calls vhags. jdvw must be at least nlon.
        !
        !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
        !            appears in the program that calls vhags. mdab must be at
        !            least min0(nlat, nlon/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
        !            appears in the program that calls vhags. ndab must be at
        !            least nlat.
        !
        !     wvhags an array which must be initialized by subroutine vhgsi.
        !            once initialized, wvhags can be used repeatedly by vhags
        !            as long as nlon and nlat remain unchanged.  wvhags must
        !            not be altered between calls of vhags.
        !
        !     lvhags the dimension of the array wvhags as it appears in the
        !            program that calls vhags. define
        !
        !               l1 = min0(nlat, nlon/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lvhags must be at least
        !
        !            (nlat+1)*(nlat+1)*nlat/2 + nlon + 15
        !
        !
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls vhags. define
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            if ityp .le. 2 then lwork must be at least
        !            the larger of the two quantities
        !
        !               3*nlat*(nlat+1)+2  (required by vhagsi)
        !
        !            and
        !
        !               (2*nt+1)*nlat*nlon
        !
        !            if ityp .gt. 2 then lwork must be at least
        !            the larger of the two quantities
        !
        !               3*nlat*(nlat+1)+2  (required by vhagsi)
        !
        !            and
        !
        !              (2*nt+1)*l2*nlon
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     br, bi  two or three dimensional arrays (see input parameter nt)
        !     cr, ci  that contain the vector spherical harmonic coefficients
        !            in the spectral representation of v(i, j) and w(i, j) given
        !            in the discription of subroutine vhses. br(mp1, np1),
        !            bi(mp1, np1), cr(mp1, np1), and ci(mp1, np1) are computed
        !            for mp1=1, ..., mmax and np1=mp1, ..., nlat except for np1=nlat
        !            and odd mp1. mmax=min0(nlat, nlon/2) if nlon is even or
        !            mmax=min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of ityp
        !            = 4  error in the specification of nt
        !            = 5  error in the specification of idvw
        !            = 6  error in the specification of jdvw
        !            = 7  error in the specification of mdab
        !            = 8  error in the specification of ndab
        !            = 9  error in the specification of lvhags
        !            = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)  :: this
        real (WP),        intent (in)      :: vector_function(:, :, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)           :: error_flag
        real (WP), allocatable :: polar_component(:, :)
        real (WP), allocatable :: azimuthal_component(:, :)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Allocate arrays
            allocate ( &
                polar_component(     1:nlat, 1:nlon), &
                azimuthal_component( 1:nlat, 1:nlon), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (sphere_t)'
                write( stderr, '(A)' ) 'Allocation failed in PERFORM_VECTOR_ANALYSIS'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! Compute the spherical angle components
        !--------------------------------------------------------------------------------

        associate( &
            F       => vector_function, &
            F_theta => polar_component, &
            F_phi   => azimuthal_component &
            )

            call this%Get_spherical_angle_components( F, F_theta, F_phi )

        end associate

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            ityp   => this%VECTOR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            v      => polar_component, &
            w      => azimuthal_component, &
            idvw   => size( polar_component, dim = 1 ), &
            jdvw   => size( polar_component, dim = 2 ), &
            br     => this%workspace%real_polar_harmonic_coefficients, &
            bi     => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr     => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci     => this%workspace%imaginary_azimuthal_harmonic_coefficients, &
            mdab   => size( this%workspace%real_polar_harmonic_coefficients, dim = 1 ), &
            ndab   => size( this%workspace%real_polar_harmonic_coefficients, dim = 2 ), &
            wvhags => this%workspace%wvhags, &
            lvhags => size( this%workspace%wvhags ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            call Vhags( nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                mdab, ndab, wvhags, lvhags, work, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'ERROR: GET_VORTICITY'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of VECTOR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_SYNTHESES'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
                write( stderr, '(A)') 'size( THETA, dim = 1 )'
                write( stderr, '(A)') 'size( PHI,   dim = 1 )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
                write( stderr, '(A)') 'size( THETA, dim = 2 )'
                write( stderr, '(A)') 'size( PHI,   dim = 2 )'


            else if ( error_flag == 7 ) then

                write( stderr, '(A)') 'Invalid extent for BR or CR'
                write( stderr, '(A)') 'size( BR, dim = 1)'
                write( stderr, '(A)') 'size( CR, dim = 1)'

            else if (error_flag == 8) then

                write( stderr, '(A)') 'Invalid extent for BI or CI'
                write( stderr, '(A)') 'size( BI, dim = 2)'
                write( stderr, '(A)') 'size( CI, dim = 2)'

            else if (error_flag == 9) then

                write( stderr, '(A)') 'Invalid extent for WVHAGS'
                write( stderr, '(A)') 'size( WVHAGS )'

            else if (error_flag == 10) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Deallocate arrays
        !--------------------------------------------------------------------------------

        deallocate ( &
            polar_component, &
            azimuthal_component, &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (sphere_t)'
            write( stderr, '(A)' ) 'Deallocation failed in PERFORM_VECTOR_ANALYSIS'
            write( stderr, '(A)' ) trim( error_message )

        end if


    end subroutine Perform_vector_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_vector_synthesis( this, polar_component, azimuthal_component )
        !
        ! Reference: http://www2.cisl.ucar.edu/spherepack/documentation#vhsgs.html
        !
        ! Documentation: SPHEREPACK 3.2
        !
        !     subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci,
        !                     mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror)
        !
        !
        !     subroutine vhsgs performs the vector spherical harmonic synthesis
        !     of the arrays br, bi, cr, and ci and stores the result in the
        !     arrays v and w.  the synthesis is performed on an equally spaced
        !     longitude grid and a gaussian colatitude grid (measured from
        !     the north pole). v(i, j) and w(i, j) are the colatitudinal and
        !     east longitudinal components respectively, located at the i(th)
        !     colatitude gaussian point (see nlat below) and longitude
        !     phi(j) = (j-1)*2*pi/nlon.  the spectral respresentation of (v, w)
        !     is given below at output parameters v, w.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0, pi) and are computed
        !            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
        !            if nlat is odd the equator will be included as the grid point
        !            theta((nlat+1)/2).  if nlat is even the equator will be
        !            excluded as a grid point and will lie half way between
        !            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
        !            note: on the half sphere, the number of grid points in the
        !            colatitudinal direction is nlat/2 if nlat is even or
        !            (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than zero. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !     ityp   = 0  no symmetries exist about the equator. the synthesis
        !                 is performed on the entire sphere.  i.e. on the
        !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon.
        !
        !            = 1  no symmetries exist about the equator. the synthesis
        !                 is performed on the entire sphere.  i.e. on the
        !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon. the curl of (v, w) is zero. that is,
        !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
        !                 the coefficients cr and ci are zero.
        !
        !            = 2  no symmetries exist about the equator. the synthesis
        !                 is performed on the entire sphere.  i.e. on the
        !                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e.,
        !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
        !                 the coefficients br and bi are zero.
        !
        !            = 3  v is symmetric and w is antisymmetric about the
        !                 equator. the synthesis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the synthesis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the synthesis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !            = 4  v is symmetric and w is antisymmetric about the
        !                 equator. the synthesis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the synthesis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the synthesis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the curl of (v, w) is zero. that is,
        !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
        !                 the coefficients cr and ci are zero.
        !
        !            = 5  v is symmetric and w is antisymmetric about the
        !                 equator. the synthesis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the synthesis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the synthesis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the divergence of (v, w) is zero. i.e.,
        !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
        !                 the coefficients br and bi are zero.
        !
        !            = 6  v is antisymmetric and w is symmetric about the
        !                 equator. the synthesis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the synthesis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the synthesis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !
        !            = 7  v is antisymmetric and w is symmetric about the
        !                 equator. the synthesis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the synthesis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the synthesis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the curl of (v, w) is zero. that is,
        !                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
        !                 the coefficients cr and ci are zero.
        !
        !            = 8  v is antisymmetric and w is symmetric about the
        !                 equator. the synthesis is performed on the northern
        !                 hemisphere only.  i.e., if nlat is odd the synthesis
        !                 is performed on the arrays v(i, j), w(i, j) for
        !                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
        !                 even the synthesis is performed on the the arrays
        !                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
        !                 the divergence of (v, w) is zero. i.e.,
        !                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
        !                 the coefficients br and bi are zero.
        !
        !
        !     nt     the number of syntheses.  in the program that calls vhsgs,
        !            the arrays v, w, br, bi, cr, and ci can be three dimensional
        !            in which case multiple syntheses will be performed.
        !            the third index is the synthesis index which assumes the
        !            values k=1, ..., nt.  for a single synthesis set nt=1. the
        !            discription of the remaining parameters is simplified
        !            by assuming that nt=1 or that all the arrays are two
        !            dimensional.
        !
        !     idvw   the first dimension of the arrays v, w as it appears in
        !            the program that calls vhags. if ityp .le. 2 then idvw
        !            must be at least nlat.  if ityp .gt. 2 and nlat is
        !            even then idvw must be at least nlat/2. if ityp .gt. 2
        !            and nlat is odd then idvw must be at least (nlat+1)/2.
        !
        !     jdvw   the second dimension of the arrays v, w as it appears in
        !            the program that calls vhsgs. jdvw must be at least nlon.
        !
        !     br, bi  two or three dimensional arrays (see input parameter nt)
        !     cr, ci  that contain the vector spherical harmonic coefficients
        !            in the spectral representation of v(i, j) and w(i, j) given
        !            below at the discription of output parameters v and w.
        !
        !     mdab   the first dimension of the arrays br, bi, cr, and ci as it
        !            appears in the program that calls vhsgs. mdab must be at
        !            least min0(nlat, nlon/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !     ndab   the second dimension of the arrays br, bi, cr, and ci as it
        !            appears in the program that calls vhsgs. ndab must be at
        !            least nlat.
        !
        !     wvhsgs an array which must be initialized by subroutine vhsgsi.
        !            once initialized, wvhsgs can be used repeatedly by vhsgs
        !            as long as nlon and nlat remain unchanged.  wvhsgs must
        !            not be altered between calls of vhsgs.
        !
        !     lvhsgs the dimension of the array wvhsgs as it appears in the
        !            program that calls vhsgs. define
        !
        !               l1 = min0(nlat, nlon/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lvhsgs must be at least
        !
        !                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat
        !
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls vhsgs. define
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            if ityp .le. 2 then lwork must be at least
        !
        !                       (2*nt+1)*nlat*nlon
        !
        !            if ityp .gt. 2 then lwork must be at least
        !
        !                        (2*nt+1)*l2*nlon
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     v, w    two or three dimensional arrays (see input parameter nt)
        !            in which the synthesis is stored. v is the colatitudinal
        !            component and w is the east longitudinal component.
        !            v(i, j), w(i, j) contain the components at the guassian colatitude
        !            point theta(i) and longitude phi(j) = (j-1)*2*pi/nlon.
        !            the index ranges are defined above at the input parameter
        !            ityp. v and w are computed from the formulas given below.
        !
        !
        !     define
        !
        !     1.  theta is colatitude and phi is east longitude
        !
        !     2.  the normalized associated legendre funnctions
        !
        !         pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
        !                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
        !                        factorial(n)) times the (n+m)th derivative
        !                        of (x**2-1)**n with respect to x=cos(theta)
        !
        !     3.  vbar(m, n, theta) = the derivative of pbar(m, n, theta) with
        !                           respect to theta divided by the square
        !                           root of n(n+1).
        !
        !         vbar(m, n, theta) is more easily computed in the form
        !
        !         vbar(m, n, theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1, n, theta)
        !         -sqrt((n-m)*(n+m+1))*pbar(m+1, n, theta))/(2*sqrt(n*(n+1)))
        !
        !     4.  wbar(m, n, theta) = m/(sin(theta))*pbar(m, n, theta) divided
        !                           by the square root of n(n+1).
        !
        !         wbar(m, n, theta) is more easily computed in the form
        !
        !         wbar(m, n, theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1))
        !         *pbar(m-1, n-1, theta)+sqrt((n-m)*(n-m-1))*pbar(m+1, n-1, theta))
        !         /(2*sqrt(n*(n+1)))
        !
        !
        !    the colatitudnal dependence of the normalized surface vector
        !                spherical harmonics are defined by
        !
        !     5.    bbar(m, n, theta) = (vbar(m, n, theta), i*wbar(m, n, theta))
        !
        !     6.    cbar(m, n, theta) = (i*wbar(m, n, theta), -vbar(m, n, theta))
        !
        !
        !    the coordinate to index mappings
        !
        !     7.   theta(i) = i(th) gaussian grid point and phi(j) = (j-1)*2*pi/nlon
        !
        !
        !     the maximum (plus one) longitudinal wave number
        !
        !     8.     mmax = min0(nlat, nlon/2) if nlon is even or
        !            mmax = min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !    if we further define the output vector as
        !
        !     9.    h(i, j) = (v(i, j), w(i, j))
        !
        !    and the complex coefficients
        !
        !     10.   b(m, n) = cmplx(br(m+1, n+1), bi(m+1, n+1))
        !
        !     11.   c(m, n) = cmplx(cr(m+1, n+1), ci(m+1, n+1))
        !
        !
        !    then for i=1, ..., nlat and  j=1, ..., nlon
        !
        !        the expansion for real h(i, j) takes the form
        !
        !     h(i, j) = the sum from n=1 to n=nlat-1 of the real part of
        !
        !         .5*(b(0, n)*bbar(0, n, theta(i))+c(0, n)*cbar(0, n, theta(i)))
        !
        !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
        !     n=nlat-1 of the real part of
        !
        !              b(m, n)*bbar(m, n, theta(i))*exp(i*m*phi(j))
        !             +c(m, n)*cbar(m, n, theta(i))*exp(i*m*phi(j))
        !
        !   *************************************************************
        !
        !   in terms of real variables this expansion takes the form
        !
        !             for i=1, ..., nlat and  j=1, ..., nlon
        !
        !     v(i, j) = the sum from n=1 to n=nlat-1 of
        !
        !               .5*br(1, n+1)*vbar(0, n, theta(i))
        !
        !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
        !     n=nlat-1 of the real part of
        !
        !       (br(m+1, n+1)*vbar(m, n, theta(i))-ci(m+1, n+1)*wbar(m, n, theta(i)))
        !                                          *cos(m*phi(j))
        !      -(bi(m+1, n+1)*vbar(m, n, theta(i))+cr(m+1, n+1)*wbar(m, n, theta(i)))
        !                                          *sin(m*phi(j))
        !
        !    and for i=1, ..., nlat and  j=1, ..., nlon
        !
        !     w(i, j) = the sum from n=1 to n=nlat-1 of
        !
        !              -.5*cr(1, n+1)*vbar(0, n, theta(i))
        !
        !     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
        !     n=nlat-1 of the real part of
        !
        !      -(cr(m+1, n+1)*vbar(m, n, theta(i))+bi(m+1, n+1)*wbar(m, n, theta(i)))
        !                                          *cos(m*phi(j))
        !      +(ci(m+1, n+1)*vbar(m, n, theta(i))-br(m+1, n+1)*wbar(m, n, theta(i)))
        !                                          *sin(m*phi(j))
        !
        !
        !      br(m+1, nlat), bi(m+1, nlat), cr(m+1, nlat), and ci(m+1, nlat) are
        !      assumed zero for m even.
        !
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of ityp
        !            = 4  error in the specification of nt
        !            = 5  error in the specification of idvw
        !            = 6  error in the specification of jdvw
        !            = 7  error in the specification of mdab
        !            = 8  error in the specification of ndab
        !            = 9  error in the specification of lvhsgs
        !            = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)   :: this
        real (WP),        intent (out)      :: polar_component(:, :)
        real (WP),        intent (out)      :: azimuthal_component(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            ityp   => this%VECTOR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            v      => polar_component, &
            w      => azimuthal_component, &
            idvw   => size( polar_component, dim = 1 ),  &
            jdvw   => size( polar_component, dim = 2 ),  &
            br     => this%workspace%real_polar_harmonic_coefficients, &
            bi     => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr     => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci     => this%workspace%imaginary_azimuthal_harmonic_coefficients, &
            mdab   => size( this%workspace%real_polar_harmonic_coefficients, dim = 1 ), &
            ndab   => size( this%workspace%real_polar_harmonic_coefficients, dim = 2 ), &
            wvhsgs => this%workspace%wvhsgs, &
            lvhsgs => size( this%workspace%wvhsgs ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            call Vhsgs( nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'ERROR: GET_VORTICITY'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of VECTOR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_SYNTHESES'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
                write( stderr, '(A)') 'size( THETA, dim = 1 )'
                write( stderr, '(A)') 'size( PHI,   dim = 1 )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'POLAR_COMPONENT (THETA) or AZIMUTHAL_COMPONENT (PHI)'
                write( stderr, '(A)') 'size( THETA, dim = 2 )'
                write( stderr, '(A)') 'size( PHI,   dim = 2 )'


            else if ( error_flag == 7 ) then

                write( stderr, '(A)') 'Invalid extent for BR or CR'
                write( stderr, '(A)') 'size( BR, dim = 1)'
                write( stderr, '(A)') 'size( CR, dim = 1)'

            else if (error_flag == 8) then

                write( stderr, '(A)') 'Invalid extent for BI or CI'
                write( stderr, '(A)') 'size( BI, dim = 2)'
                write( stderr, '(A)') 'size( CI, dim = 2)'

            else if (error_flag == 9) then

                write( stderr, '(A)') 'Invalid extent for WVHSGS'
                write( stderr, '(A)') 'size( WVHSGS )'

            else if (error_flag == 10) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Perform_vector_synthesis
    !
    !*****************************************************************************************
    ! TODO
    subroutine Get_legendre_functions( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Get_legendre_functions
    !
    !*****************************************************************************************
    ! TODO
    subroutine Icosahedral_geodesic( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Icosahedral_geodesic
    !
    !*****************************************************************************************
    ! TODO
    subroutine Perform_multiple_ffts( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Perform_multiple_ffts
    !
    !*****************************************************************************************
    !
    subroutine Get_gaussian_weights_and_points( this, nlat, theta, wts )
        !
        !< Purpose:
        !
        ! Computes the nlat-many gaussian (co)latitudes and weights.
        ! the colatitudes are in radians and lie in the interval (0, pi).
        !
        ! References:
        !
        ! [1] Swarztrauber, Paul N.
        !     "On computing the points and weights for Gauss--Legendre quadrature."
        !     SIAM Journal on Scientific Computing 24.3 (2003): 945-954.

        ! [2]  http://www2.cisl.ucar.edu/spherepack/documentation#gaqd.html
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t),       intent (in out)   :: this
        integer (IP),           intent (in)       :: nlat     !! number of latitudinal points
        real (WP), allocatable, intent (out)      :: theta(:) !! latitudinal points: 0 <= theta <= pi
        real (WP), allocatable, intent (out)      :: wts(:)   !! gaussian weights
        !--------------------------------------------------------------------------------

        call this%grid%Get_gaussian_weights_and_points( nlat, theta, wts )

    end subroutine Get_gaussian_weights_and_points
    !
    !*****************************************************************************************
    !
    subroutine Finalize_sphere( this )
        !
        !< Purpose:
        !< Finalize object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%Destroy()

    end subroutine Finalize_sphere
    !
    !*****************************************************************************************
    !
end module type_sphere_mod
