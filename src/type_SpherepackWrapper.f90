module type_SpherepackWrapper

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    use type_SpherepackWorkspace, only: &
        SpherepackWorkspace

    use type_GaussianGrid, only: &
        GaussianGrid

    use type_SphericalUnitVectors, only: &
        SphericalUnitVectors

    use type_TrigonometricFunctions, only: &
        TrigonometricFunctions

    use type_ThreeDimensionalVector, only: &
        ThreeDimensionalVector, &
        assignment(=), &
        operator(*)
    
    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: SpherepackWrapper

    !----------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------------------------------------
    character (len=250) :: error_message     !! Probably long enough
    integer (ip)        :: allocate_status   !! To check allocation status
    integer (ip)        :: deallocate_status !! To check deallocation status
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, public :: SpherepackWrapper
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                       public  :: initialized = .false. !! Instantiation status
        integer (ip),                  public  :: NLON = 0   !! number of longitudinal points
        integer (ip),                  public  :: NLAT = 0   !! number of latitudinal points
        integer (ip),                  public  :: NTRUNC = 0 !! triangular truncation limit
        integer (ip),                  public  :: SCALAR_SYMMETRIES = 0 !! symmetries about the equator for scalar calculations
        integer (ip),                  public  :: VECTOR_SYMMETRIES = 0 !! symmetries about the equator for vector calculations
        integer (ip),                  public  :: NUMBER_OF_SYNTHESES = 0
        complex (wp), allocatable,     public  :: complex_spectral_coefficients(:)
        type (SpherepackWorkspace),    private :: workspace
        type (GaussianGrid),           public  :: grid
        type (TrigonometricFunctions), public  :: trigonometric_functions
        type (SphericalUnitVectors),   public  :: unit_vectors
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public  :: get_colatitude_derivative !! Vtsgs
        procedure, public  :: get_gradient !! Gradgs
        procedure, public  :: invert_gradient !!  Igradgs
        procedure, public  :: get_divergence !! Divgs
        procedure, public  :: invert_divergence !!Idivgs
        procedure, public  :: get_vorticity !! Vrtgs
        procedure, public  :: invert_vorticity !! Ivrtgs
        procedure, public  :: invert_divergence_and_vorticity !! Idvtgs
        procedure, public  :: get_scalar_laplacian !! Slapgs
        procedure, public  :: invert_helmholtz !! Islapgs
        procedure, public  :: get_vector_laplacian !! Vlapgs
        procedure, public  :: invert_vector_laplacian !! Ivlapgs
        procedure, public  :: get_stream_function_and_velocity_potential
        procedure, public  :: invert_stream_function_and_velocity_potential
        procedure, public  :: perform_grid_transfers
        procedure, public  :: perform_geo_math_coordinate_transfers
        procedure, public  :: perform_scalar_analysis
        procedure, public  :: perform_scalar_synthesis
        procedure, public  :: perform_scalar_projection !! Shpg
        procedure, public  :: perform_vector_analysis
        procedure, public  :: perform_vector_synthesis
        procedure, public  :: get_Legendre_functions
        procedure, public  :: perform_complex_analysis
        procedure, public  :: perform_complex_synthesis
        procedure, public  :: create => create_spherepack_wrapper
        procedure, public  :: destroy => destroy_spherepack_wrapper
        procedure, public  :: create_spherepack_wrapper
        procedure, public  :: destroy_spherepack_wrapper
        procedure, public  :: get_index
        procedure, public  :: get_coefficient
        procedure, public  :: compute_surface_integral
        procedure, public  :: get_rotation_operator
        procedure, public  :: synthesize_from_spec
        procedure, private :: assert_initialized
        procedure, private :: get_scalar_symmetries
        procedure, private :: get_vector_symmetries
        final              :: finalize_spherepack_wrapper
        !----------------------------------------------------------------------
    end type SpherepackWrapper


contains


    subroutine create_spherepack_wrapper( this, nlat, nlon, isym, itype )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)  :: this
        integer (ip),     intent (in)               :: nlat
        integer (ip),     intent (in)               :: nlon
        integer (ip),     intent (in), optional    :: isym      !! Either 0, 1, or 2
        integer (ip),     intent (in), optional    :: itype     !! Either 0, 1, 2, 3, ..., 8
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy_spherepack_wrapper()

        ! Set constants
        this%NLAT = nlat
        this%NLON = nlon
        this%NTRUNC = nlat - 1 !! Set triangular truncation
        this%NUMBER_OF_SYNTHESES = 1

        ! Set scalar symmetries
        if ( present( isym ) ) then
            call this%get_scalar_symmetries( isym )
        end if

        ! Set vector symmetries
        if (present (itype) ) then
            call this%get_vector_symmetries( itype )
        end if

        !--------------------------------------------------------------------------------
        ! Allocate array
        !--------------------------------------------------------------------------------

        associate( size_spec => nlat * (nlat + 1)/2 )

            ! Allocate pointer for complex spectral coefficients
            allocate( &
                this%complex_spectral_coefficients( 1:size_spec ), &
                stat=allocate_status, &
                errmsg = error_message )

            ! Check allocate status
            if ( allocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
                write( stderr, '(A)' ) 'Allocating COMPLEX_SPECTRAL_COEFFICIENTS failed in CREATE_spherepack_wrapper'
                write( stderr, '(A)' ) trim( error_message )
            end if

        end associate

        !--------------------------------------------------------------------------------
        ! create derived data types
        !--------------------------------------------------------------------------------

        call this%grid%create( nlat, nlon )

        call this%workspace%create( nlat, nlon )

        ! Set grids to compute frequently used trigonometric functions
        associate( &
            theta => this%grid%latitudes, &
            phi => this%grid%longitudes &
            )

            call this%trigonometric_functions%create( theta, phi )

        end associate

        ! Compute spherical unit vectors
        associate( &
            sint => this%trigonometric_functions%sint, &
            cost => this%trigonometric_functions%cost, &
            sinp => this%trigonometric_functions%sinp, &
            cosp => this%trigonometric_functions%cosp &
            )

            call this%unit_vectors%create( sint, cost, sinp, cosp )

        end associate

        ! Set initialization flag
        this%initialized = .true.
        
    end subroutine create_spherepack_wrapper
    !
    
    !
    subroutine destroy_spherepack_wrapper( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check status
        !--------------------------------------------------------------------------------

        if ( .not. this%initialized ) return

        !--------------------------------------------------------------------------------
        ! Deallocate array
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%complex_spectral_coefficients ) ) then

            ! Deallocate array
            deallocate( &
                this%complex_spectral_coefficients, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
                write( stderr, '(A)' ) 'Deallocating COMPLEX_SPECTRAL_COEFFICIENTS failed in DESTROY_spherepack_wrapper'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        !--------------------------------------------------------------------------------
        ! destroy derived data types
        !--------------------------------------------------------------------------------

        call this%grid%destroy()
        call this%workspace%destroy()
        call this%trigonometric_functions%destroy()
        call this%unit_vectors%destroy()

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
        ! Reset initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .false.

    end subroutine destroy_spherepack_wrapper
    !
    
    !
    subroutine assert_initialized( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. this%initialized ) then

            write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
            write( stderr, '(A)' ) 'You must instantiate object before calling methods'

        end if

    end subroutine assert_initialized
    !
    
    !
    function get_index( this, n, m ) result( return_value )
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
        ! integer (ip), dimension ((MTRUNC+1)*(MTRUNC+2)/2) :: indxm, indxn
        ! indxm = [((m, n=m, MTRUNC), m=0, MTRUNC)]
        ! indxn = [((n, n=m, MTRUNC), m=0, MTRUNC)]
        !
        ! Conversely, the index nm as a function of m and n is:
        ! nm = sum([(i, i=MTRUNC+1, MTRUNC-m+2, -1)])+n-m+1
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)  :: this
        integer (ip),     intent (in)      :: n
        integer (ip),     intent (in)      :: m
        integer (ip):: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: i !! Counter
        !--------------------------------------------------------------------------------

        associate( ntrunc => this%NTRUNC )

            if ( m <= n .and. max( n, m ) <= ntrunc ) then

                return_value = &
                    sum ( [ (i, i = ntrunc+1, ntrunc - m + 2, - 1) ] ) + n - m + 1
            else

                return_value = -1

            end if

        end associate

    end function get_index
    !
    
    !
    function get_coefficient( this, n, m ) result ( return_value )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)    :: this
        integer (ip),     intent (in)        :: n
        integer (ip),     intent (in)        :: m
        complex (wp):: return_value
        !--------------------------------------------------------------------------------

        associate( &
            ntrunc   => this%NTRUNC, &
            nm       => this%get_index( n, m ), &
            nm_conjg => this%get_index(n, -m ), &
            psi      => this%complex_spectral_coefficients &
            )

            if ( m < 0 .and. nm_conjg > 0 ) then
                return_value = ( (-1.0_wp)**(-m) ) * conjg( psi(nm_conjg) )
            else if ( nm > 0 ) then
                return_value = psi(nm)
            else
                return_value = 0.0_wp
            end if

        end associate

    end function get_coefficient
    !
    
    !
    subroutine get_scalar_symmetries( this, isym )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        integer (ip),              intent (in)     :: isym
        !--------------------------------------------------------------------------------

        select case (isym)
            case (2)
        		
                this%SCALAR_SYMMETRIES = isym
            case (1)
        		
                this%SCALAR_SYMMETRIES = isym
            case (0)
        		
                this%SCALAR_SYMMETRIES = isym
            case default
        		
                ! Handle invalid symmetry arguments
                write( stderr, '(A)' )     'TYPE (SpherepackWrapper) in get_SCALAR_SYMMETRIES'
                write( stderr, '(A, I2)' ) 'Optional argument isym = ', isym
                write( stderr, '(A)' )     'must be either 0, 1, or 2 (default isym = 0)'
        end select

    end subroutine get_scalar_symmetries
        !
    
    !
    subroutine get_vector_symmetries( this, itype )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)  :: this
        integer (ip),     intent (in)      :: itype
        !--------------------------------------------------------------------------------

        select case (itype)
            case (8)
                this%VECTOR_SYMMETRIES = itype
            case (7)
                this%VECTOR_SYMMETRIES = itype
            case (6)
                this%VECTOR_SYMMETRIES = itype
            case (5)
                this%VECTOR_SYMMETRIES = itype
            case (4)
                this%VECTOR_SYMMETRIES = itype
            case (3)
                this%VECTOR_SYMMETRIES = itype
            case (2)
                this%VECTOR_SYMMETRIES = itype
            case (1)
                this%VECTOR_SYMMETRIES = itype
            case (0)
                this%VECTOR_SYMMETRIES = itype
            case default
                ! Handle invalid symmetry arguments
                write( stderr, '(A)' )     'TYPE (SpherepackWrapper) in get_VECTOR_SYMMETRIES'
                write( stderr, '(A, I2)' ) 'Optional argument itype = ', itype
                write( stderr, '(A)' )     'must be either 0, 1, 2, ..., 8 (default itype = 0)'
        end select

    end subroutine get_vector_symmetries
    !
    
    !
    subroutine perform_complex_analysis( this, scalar_function )
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
        ! integer (ip), dimension ((mtrunc+1)*(mtrunc+2)/2) :: indxm, indxn
        ! indxm = [((m, n=m, mtrunc), m=0, mtrunc)]
        ! indxn = [((n, n=m, mtrunc), m=0, mtrunc)]
        !
        ! Conversely, the index nm as a function of m and n is:
        ! nm = sum([(i, i=mtrunc+1, mtrunc-m+2, -1)])+n-m+1
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)  :: this
        real (wp),                 intent (in)      :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: m  !< counters
        integer (ip)::  n  !< counters
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()
        
        !--------------------------------------------------------------------------------
        ! compute the (real) spherical harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_scalar_analysis( scalar_function )

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
                0.5_wp * [((a(m + 1, n + 1), n = m, ntrunc), m = 0, ntrunc)], &
                0.5_wp * [((b(m + 1, n + 1), n = m, ntrunc), m = 0, ntrunc)], &
                wp )
 
        end associate

    end subroutine perform_complex_analysis
    !
    
    !
    subroutine perform_complex_synthesis( this, scalar_function )
        !
        !< Purpose:
        ! converts gridded input array (datagrid) to (complex) spherical harmonic coefficients
        ! (dataspec).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)  :: this
        real (wp),                 intent (out)     :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: m  !! Counters
        integer (ip):: n  !! Counters
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

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
                    associate(nm => this%get_index( n, m ))
                
                        ! set the real component
                        a( m + 1, n + 1 ) = 2.0_wp * real( psi(nm) )
                
                        ! set the imaginary component
                        b( m + 1, n + 1 ) = 2.0_wp * aimag( psi(nm) )

                    end associate
                end do
            end do
        
        end associate

        !--------------------------------------------------------------------------------
        ! synthesise the scalar function from the (real) harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_scalar_synthesis( scalar_function )
 
    end subroutine perform_complex_synthesis
    !
    
    !
    subroutine synthesize_from_spec( this, spec, scalar_function )
        !
        !< Purpose:
        !
        ! Used mainly for testing
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)  :: this
        complex (wp),              intent (in)      :: spec(:)
        real (wp),                 intent (out)     :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: m  !! Counters
        integer (ip):: n  !! Counters
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

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
                    spec_index: associate( nm => this%get_index( n, m ) )
                
                        ! set the real component
                        a(m + 1, n + 1) = 2.0_wp * real( spec(nm) )
                
                        ! set the imaginary component
                        b(m + 1, n + 1) = 2.0_wp * aimag( spec(nm) )

                    end associate spec_index

                end do
            end do
        
        end associate

        !--------------------------------------------------------------------------------
        ! synthesise the scalar function from the (real) harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_scalar_synthesis( scalar_function )
 
    end subroutine synthesize_from_spec
    !
    
    !
    subroutine get_rotation_operator( this, scalar_function, rotation_operator)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)  :: this
        real (wp),                 intent (in)      :: scalar_function(:, :)
        real (wp),                 intent (out)     :: rotation_operator(:, :, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)            :: k, l   !! Counters
        real (wp), allocatable :: polar_gradient_component(:, :)
        real (wp), allocatable :: azimuthal_gradient_component(:, :)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Allocate arrays
            allocate( &
                polar_gradient_component(     nlat, nlon ), &
                azimuthal_gradient_component( nlat, nlon ), &
                stat=allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
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

            call this%get_gradient( f, grad_theta, grad_phi )

        end associate

        !--------------------------------------------------------------------------------
        ! Calculate the rotation operator applied to a scalar function
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            do l = 1, nlon
                do k = 1, nlat

                    associate( &
                        theta      => this%unit_vectors%polar(k, l), &
                        phi        => this%unit_vectors%azimuthal( k, l), &
                        grad_theta => polar_gradient_component(k, l), &
                        grad_phi   => azimuthal_gradient_component(k, l) &
                        )

                        rotation_operator(:, k, l) = phi * grad_theta - theta * grad_phi

                    end associate
                end do
            end do
        end associate

        !--------------------------------------------------------------------------------
        ! Release memory
        !--------------------------------------------------------------------------------

        deallocate( &
            polar_gradient_component, &
            azimuthal_gradient_component, &
            stat=deallocate_status, &
            errmsg = error_message )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
            write( stderr, '(A)' ) 'Deallocation failed in GET_ROTATION_OPERATOR'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end subroutine get_rotation_operator
    !
    
    !
    function compute_surface_integral( this, scalar_function ) result( return_value )
        !
        !< Purpose:
        !
        ! computes the (scalar) surface integral on the sphere (S^2):
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
        class (SpherepackWrapper), intent (in out) :: this
        real (wp),                 intent (in)     :: scalar_function(:, :)
        real (wp)                                   :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)            :: k            !! counter
        real (wp), allocatable :: summation(:)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate array
        !--------------------------------------------------------------------------------

        associate( nlat => this%NLAT )

            ! Allocate array
            allocate( &
                summation( nlat ), &
                stat=allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
                write( stderr, '(A)' ) 'Allocation failed in COMPUTE_SURFACE_INTEGRAL'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! compute integral
        !--------------------------------------------------------------------------------

        ! Initialize array
        summation = 0.0_wp

        ! compute the integrant
        associate( &
            nlat => this%NLAT, &
            Dphi => this%grid%mesh_phi, &
            wts  => this%grid%gaussian_weights, &
            f    => scalar_function &
            )

            ! Apply trapezoidal rule
            do k = 1, nlat

                summation(k) = sum( f(k, :) ) * dphi

            end do

            ! Apply gaussian quadrature
            summation = summation * wts

        end associate

        ! Set integral \int_{S^2} f( theta, phi ) dS
        return_value = sum( summation )

        !--------------------------------------------------------------------------------
        ! Deallocate array
        !--------------------------------------------------------------------------------

        deallocate( &
            summation, &
            stat=deallocate_status, &
            errmsg = error_message )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then
            write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
            write( stderr, '(A)' ) 'Deallocation failed in COMPUTE_SURFACE_INTEGRAL'
            write( stderr, '(A)' ) trim( error_message )
        end if

    end function compute_surface_integral
    !
    
    !
    subroutine compute_first_moment( this, scalar_function, first_moment )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper),      intent (in out)  :: this
        real (wp),                      intent (in)      :: scalar_function(:, :)
        type (ThreeDimensionalVector),  intent (out)     :: first_moment
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)            :: k, l !! Counters
        real (wp), allocatable :: integrant(:, :, :)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate array
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Allocate arrays
            allocate( &
                integrant( nlat, nlon, 3 ), &
                stat=allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
                write( stderr, '(A)' ) 'Allocation failed in COMPUTE_FIRST_MOMENT'
                write( stderr, '(A)' ) trim( error_message )
            end if

        end associate

        !--------------------------------------------------------------------------------
        ! compute integrant
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! compute integrant
            do l = 1, nlon
                do k = 1, nlat

                    associate( &
                        u => this%unit_vectors%radial( k, l ), &
                        f => scalar_function(k, l) &
                        )

                        integrant(k, l, 1) = u%x * f
                        integrant(k, l, 2) = u%y * f
                        integrant(k, l, 3) = u%z * f

                    end associate
                end do
            end do
        end associate

        !--------------------------------------------------------------------------------
        ! compute first moment
        !--------------------------------------------------------------------------------


        associate( &
            M  => first_moment, &
            f1 => integrant(:, :, 1), &
            f2 => integrant(:, :, 2), &
            f3 => integrant(:, :, 3) &
            )

            m%x = this%compute_surface_integral( f1 )

            m%y = this%compute_surface_integral( f2 )

            m%z = this%compute_surface_integral( f3 )

        end associate

        !--------------------------------------------------------------------------------
        ! Deallocate array
        !--------------------------------------------------------------------------------

        deallocate( &
            integrant, &
            stat=deallocate_status, &
            errmsg = error_message )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
            write( stderr, '(A)' ) 'Deallocation failed in COMPUTE_FIRST_MOMENT'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end subroutine compute_first_moment
    !
    
    !
    ! Public SPHEREPACK 3.2 methods
    !
    
    !
    subroutine get_colatitude_derivative( this, polar_component, azimuthal_component )
        !
        !< Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vtsgs.html
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)   :: this
        real (wp),        intent (out)      :: polar_component(:)     !! vt
        real (wp),        intent (out)      :: azimuthal_component(:) !! wt
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

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

    end subroutine get_colatitude_derivative
    !
    
    !
    subroutine get_gradient( this, scalar_function, &
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
        class (SpherepackWrapper), intent (in out)  :: this
        real (wp),        intent (in)      :: scalar_function(:, :)
        real (wp),        intent (out)     :: polar_gradient_component(:, :)
        real (wp),        intent (out)     :: azimuthal_gradient_component(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! compute the (real) spherical harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_scalar_analysis( scalar_function )

        !--------------------------------------------------------------------------------
        ! compute gradient
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

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: get_GRADIENT'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of SCALAR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'

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

    end subroutine get_gradient
    !
    
    !
    subroutine invert_gradient( this, &
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
        class (SpherepackWrapper), intent (in out) :: this
        real (wp),        intent (in)     :: polar_gradient_component(:, :)
        real (wp),        intent (in)     :: azimuthal_gradient_component(:, :)
        real (wp),        intent (out)    :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

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

    end subroutine invert_gradient
    !
    
    !
    subroutine get_divergence( this, vector_field, divergence )
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
        class (SpherepackWrapper), intent (in out) :: this
        real (wp),        intent (in)     :: vector_field (:, :, :)
        real (wp),        intent (out)    :: divergence (:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Calculate the (real) vector harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_vector_analysis( vector_field )

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

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: get_DIVERGENCE'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of SCALAR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'

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


    end subroutine get_divergence
    !
    
    !
    subroutine invert_divergence( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine invert_divergence
    !
    
    !
    subroutine get_vorticity( this, vector_function, vorticity )
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
        class (SpherepackWrapper), intent (in out)  :: this
        real (wp),        intent (in)      :: vector_function(:, :, :)
        real (wp),        intent (out)     :: vorticity(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! calculate the (real) vector harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_vector_analysis( vector_function )

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

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'ERROR: get_VORTICITY'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'


            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of SCALAR_SYMMETRIES'


            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'



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

    end subroutine get_vorticity
    !
    
    !
    subroutine invert_vorticity( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine invert_vorticity
    !
    
    !
    subroutine invert_divergence_and_vorticity( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine invert_divergence_and_vorticity
    !
    
    !
    subroutine get_scalar_laplacian( this, scalar_function, scalar_laplacian )
        !
        ! References:
        ! http://www2.cisl.ucar.edu/spherepack/documentation#slapgs.html
        !
        ! Documentation: SPHEREPACK 3.2
        !
        !
        !     subroutine slapgs(nlat, nlon, isym, nt, slap, ids, jds, a, b,
        !    mdab, ndab, wshsgs, lshsgs, work, lwork, ierror)
        !
        !
        !     given the scalar spherical harmonic coefficients a and b, precomputed
        !     by subroutine shags for a scalar field sf, subroutine slapgs computes
        !     the laplacian of sf in the scalar array slap.  slap(i, j) is the
        !     laplacian of sf at the gaussian colatitude theta(i) (see nlat as
        !     an input parameter) and east longitude lambda(j) = (j-1)*2*pi/nlon
        !     on the sphere.  i.e.
        !
        !         slap(i, j) =
        !
        !                  2                2
        !         [1/sint*d (sf(i, j)/dlambda + d(sint*d(sf(i, j))/dtheta)/dtheta]/sint
        !
        !
        !     where sint = sin(theta(i)).  the scalar laplacian in slap has the
        !     same symmetry or absence of symmetry about the equator as the scalar
        !     field sf.  the input parameters isym, nt, mdab, ndab must have the
        !     same values used by shags to compute a and b for sf. the associated
        !     legendre functions are stored rather than recomputed as they are
        !     in subroutine slapgc.
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
        !     nlon   the number of distinct longitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than zero. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !     isym   this parameter should have the same value input to subroutine
        !            shags to compute the coefficients a and b for the scalar field
        !            sf.  isym is set as follows:
        !
        !            = 0  no symmetries exist in sf about the equator. scalar
        !                 synthesis is used to compute slap on the entire sphere.
        !                 i.e., in the array slap(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon.
        !
        !           = 1  sf and slap are antisymmetric about the equator. the
        !                synthesis used to compute slap is performed on the
        !                northern hemisphere only.  if nlat is odd, slap(i, j) is
        !                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
        !                nlat is even, slap(i, j) is computed for i=1, ..., nlat/2
        !                and j=1, ..., nlon.
        !
        !
        !           = 2  sf and slap are symmetric about the equator. the
        !                synthesis used to compute slap is performed on the
        !                northern hemisphere only.  if nlat is odd, slap(i, j) is
        !                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
        !                nlat is even, slap(i, j) is computed for i=1, ..., nlat/2
        !                and j=1, ..., nlon.
        !
        !
        !     nt     the number of analyses.  in the program that calls slapgs
        !            the arrays slap, a, and b can be three dimensional in which
        !            case multiple synthesis will be performed.  the third index
        !            is the synthesis index which assumes the values k=1, ..., nt.
        !            for a single analysis set nt=1. the description of the
        !            remaining parameters is simplified by assuming that nt=1
        !            or that all the arrays are two dimensional.
        !
        !   ids      the first dimension of the array slap as it appears in the
        !            program that calls slapgs.  if isym = 0 then ids must be at
        !            least nlat.  if isym > 0 and nlat is even then ids must be
        !            at least nlat/2. if isym > 0 and nlat is odd then ids must
        !            be at least (nlat+1)/2.
        !
        !   jds      the second dimension of the array slap as it appears in the
        !            program that calls slapgs. jds must be at least nlon.
        !
        !
        !   a, b      two or three dimensional arrays (see input parameter nt)
        !            that contain scalar spherical harmonic coefficients
        !            of the scalar field sf as computed by subroutine shags.
        !     ***    a, b must be computed by shags prior to calling slapgs.
        !
        !
        !    mdab    the first dimension of the arrays a and b as it appears
        !            in the program that calls slapgs.  mdab must be at
        !            least min0(nlat, (nlon+2)/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !    ndab    the second dimension of the arrays a and b as it appears
        !            in the program that calls slapgs. ndbc must be at least
        !            least nlat.
        !
        !            mdab, ndab should have the same values input to shags to
        !            compute the coefficients a and b.
        !
        !
        !    wshsgs  an array which must be initialized by subroutine slapgsi
        !            (or equivalently by shsgsi).  once initialized, wshsgs
        !            can be used repeatedly by slapgs as long as nlat and nlon
        !            remain unchanged.  wshsgs must not be altered between calls
        !            of slapgs.
        !
        !    lshsgs  the dimension of the array wshsgs as it appears in the
        !            program that calls slapgs.  let
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
        !            program that calls slapgs. define
        !
        !               l2 = nlat/2                    if nlat is even or
        !               l2 = (nlat+1)/2                if nlat is odd
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            if isym is zero then lwork must be at least
        !
        !               (nt+1)*nlat*nlon + nlat*(2*nt*l1+1)
        !
        !            if isym is nonzero lwork must be at least
        !
        !               (nt+1)*l2*nlon + nlat*(2*nt*l1+1)
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !
        !    slap    a two or three dimensional arrays (see input parameter nt) that
        !            contain the scalar laplacian of the scalar field sf.  slap(i, j)
        !            is the scalar laplacian at the gaussian colatitude theta(i)
        !            and longitude lambda(j) = (j-1)*2*pi/nlon for i=1, ..., nlat
        !            and j=1, ..., nlon.
        !
        !
        !  ierror    a parameter which flags errors in input parameters as follows:
        !
        !            = 0  no errors detected
        !
        !            = 1  error in the specification of nlat
        !
        !            = 2  error in the specification of nlon
        !
        !            = 3  error in the specification of ityp
        !
        !            = 4  error in the specification of nt
        !
        !            = 5  error in the specification of ids
        !
        !            = 6  error in the specification of jds
        !
        !            = 7  error in the specification of mdbc
        !
        !            = 8  error in the specification of ndbc
        !
        !            = 9  error in the specification of lshsgs
        !
        !            = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        real (wp),        intent (in)     :: scalar_function(:,:)
        real (wp),        intent (out)    :: scalar_laplacian(:,:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Set (real) scalar spherica harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_scalar_analysis( scalar_function )

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            isym   => this%SCALAR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            slap   => scalar_laplacian, &
            ids    => size( scalar_laplacian, dim = 1), &
            jds    => size( scalar_laplacian, dim = 2), &
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients, &
            mdab   => size( this%workspace%real_harmonic_coefficients, dim = 1 ), &
            ndab   => size( this%workspace%real_harmonic_coefficients, dim = 2 ), &
            wshsgs => this%workspace%wshsgs, &
            lshsgs => size( this%workspace%wshsgs ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            ierror => error_flag &
            )

            call Slapgs( nlat, nlon, isym, nt, slap, ids, jds, a, b, &
                mdab, ndab, wshsgs, lshsgs, work, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: get_SCALAR_LAPLACIAN'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of VECTOR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for SCALAR_LAPLACIAN'
                write( stderr, '(A)') 'size( SCALAR_LAPLACIAN, dim = 1 )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Invalid extent for SCALAR_LAPLACIAN'
                write( stderr, '(A)') 'size( SCALAR_LAPLACIAN, dim = 2 )'

            else if ( error_flag == 7 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A) '
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
                write( stderr, '(A)') 'size( SOLUTION, dim = 1 )'

            else if ( error_flag == 8 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A) '
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
                write( stderr, '(A)') 'size( SOLUTION, dim = 2 )'

            else if ( error_flag == 9 ) then

                write( stderr, '(A)') 'Invalid extent for WSHSGS'
                write( stderr, '(A)') 'size( WSHSGS )'

            else if ( error_flag == 10 ) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine get_scalar_laplacian
    !
    
    !
    subroutine invert_helmholtz( this, &
        helmholtz_constant, source_term, solution, perturbation_optional )
        !
        ! Reference:
        ! http://www2.cisl.ucar.edu/spherepack/documentation#islapgs.html
        !
        ! Documentation: SPHEREPACK 3.2
        !
        !
        !     subroutine islapgs(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b,
        !    mdab, ndab, wshsgs, lshsgs, work, lwork, pertrb, ierror)
        !
        !     islapgs inverts the laplace or helmholz operator on a Gaussian grid.
        !     Given the spherical harmonic coefficients a(m, n) and b(m, n) of the
        !     right hand side slap(i, j), islapgc computes a solution sf(i, j) to
        !     the following helmhotz equation :
        !
        !           2                2
        !     [d(sf(i, j))/dlambda /sint + d(sint*d(sf(i, j))/dtheta)/dtheta]/sint
        !
        !                   - xlmbda * sf(i, j) = slap(i, j)
        !
        !      where sf(i, j) is computed at the Gaussian colatitude point theta(i)
        !      (see nlat as an input argument) and longitude
        !
        !                 lambda(j) = (j-1)*2*pi/nlon
        !
        !            for i=1, ..., nlat and j=1, ..., nlon.
        !
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
        !     nlon   the number of distinct longitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than zero. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !     isym   this parameter should have the same value input to subroutine
        !            shags to compute the coefficients a and b for the scalar field
        !            slap.  isym is set as follows:
        !
        !            = 0  no symmetries exist in slap about the equator. scalar
        !                 synthesis is used to compute sf on the entire sphere.
        !                 i.e., in the array sf(i, j) for i=1, ..., nlat and
        !                 j=1, ..., nlon.
        !
        !           = 1  sf and slap are antisymmetric about the equator. the
        !                synthesis used to compute sf is performed on the
        !                northern hemisphere only.  if nlat is odd, sf(i, j) is
        !                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
        !                nlat is even, sf(i, j) is computed for i=1, ..., nlat/2
        !                and j=1, ..., nlon.
        !
        !
        !           = 2  sf and slap are symmetric about the equator. the
        !                synthesis used to compute sf is performed on the
        !                northern hemisphere only.  if nlat is odd, sf(i, j) is
        !                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
        !                nlat is even, sf(i, j) is computed for i=1, ..., nlat/2
        !                and j=1, ..., nlon.
        !
        !
        !     nt     the number of analyses.  in the program that calls islapgs
        !            the arrays sf, a, and b can be three dimensional in which
        !            case multiple synthesis will be performed.  the third index
        !            is the synthesis index which assumes the values k=1, ..., nt.
        !            k is also the index for the perturbation array pertrb.
        !            for a single analysis set nt=1. the description of the
        !            remaining parameters is simplified by assuming that nt=1
        !            or that sf, a, b are two dimensional and pertrb is a constant.
        !
        !   xlmbda   a one dimensional array with nt elements. if xlmbda is
        !            is identically zero islapgc solves poisson's equation.
        !            if xlmbda > 0.0 islapgc solves the helmholtz equation.
        !            if xlmbda < 0.0 the nonfatal error flag ierror=-1 is
        !            returned. negative xlambda could result in a division
        !            by zero.
        !
        !   ids      the first dimension of the array sf as it appears in the
        !            program that calls islapgs.  if isym = 0 then ids must be at
        !            least nlat.  if isym > 0 and nlat is even then ids must be
        !            at least nlat/2. if isym > 0 and nlat is odd then ids must
        !            be at least (nlat+1)/2.
        !
        !   jds      the second dimension of the array sf as it appears in the
        !            program that calls islapgs. jds must be at least nlon.
        !
        !
        !   a, b      two or three dimensional arrays (see input parameter nt)
        !            that contain scalar spherical harmonic coefficients
        !            of the scalar field slap as computed by subroutine shags.
        !     ***    a, b must be computed by shags prior to calling islapgs.
        !
        !
        !    mdab    the first dimension of the arrays a and b as it appears
        !            in the program that calls islapgs.  mdab must be at
        !            least min0(nlat, (nlon+2)/2) if nlon is even or at least
        !            min0(nlat, (nlon+1)/2) if nlon is odd.
        !
        !    ndab    the second dimension of the arrays a and b as it appears
        !            in the program that calls islapgs. ndbc must be at least
        !            least nlat.
        !
        !            mdab, ndab should have the same values input to shags to
        !            compute the coefficients a and b.
        !
        !
        !    wshsgs  an array which must be initialized by subroutine islapgsi
        !            (or equivalently by shsesi).  once initialized, wshsgs
        !            can be used repeatedly by islapgs as long as nlat and nlon
        !            remain unchanged.  wshsgs  must not be altered between calls
        !            of islapgs.
        !
        !    lshsgs  the dimension of the array wshsgs as it appears in the
        !            program that calls islapgs.  let
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
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls islapgs. define
        !
        !               l2 = nlat/2                    if nlat is even or
        !               l2 = (nlat+1)/2                if nlat is odd
        !               l1 = min0(nlat, (nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat, (nlon+1)/2) if nlon is odd
        !
        !            if isym is zero then lwork must be at least
        !
        !               (nt+1)*nlat*nlon + nlat*(2*nt*l1+1)
        !
        !            if isym is nonzero lwork must be at least
        !
        !               (nt+1)*l2*nlon + nlat*(2*nt*l1+1)
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !
        !    sf      a two or three dimensional arrays (see input parameter nt) that
        !            inverts the scalar laplacian in slap.  sf(i, j) is given at
        !            the colatitude
        !
        !                 theta(i) = (i-1)*pi/(nlat-1)
        !
        !            and longitude
        !
        !                 lambda(j) = (j-1)*2*pi/nlon
        !
        !            for i=1, ..., nlat and j=1, ..., nlon.
        !
        !   pertrb  a one dimensional array with nt elements (see input
        !           parameter nt). in the discription that follows we assume
        !           that nt=1. if xlmbda > 0.0 then pertrb=0.0 is always
        !           returned because the helmholtz operator is invertible.
        !           if xlmbda = 0.0 then a solution exists only if a(1, 1)
        !           is zero. islapec sets a(1, 1) to zero. the resulting
        !           solution sf(i, j) solves poisson's equation with
        !           pertrb = a(1, 1)/(2.*sqrt(2.)) subtracted from the
        !           right side slap(i, j).
        !
        !  ierror    a parameter which flags errors in input parameters as follows:
        !
        !            = 0  no errors detected
        !
        !            = 1  error in the specification of nlat
        !
        !            = 2  error in the specification of nlon
        !
        !            = 3  error in the specification of ityp
        !
        !            = 4  error in the specification of nt
        !
        !            = 5  error in the specification of ids
        !
        !            = 6  error in the specification of jds
        !
        !            = 7  error in the specification of mdbc
        !
        !            = 8  error in the specification of ndbc
        !
        !            = 9  error in the specification of lshsgs
        !
        !            = 10 error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)        :: this
        real (wp),        intent (in)            :: helmholtz_constant
        real (wp),        intent (in)            :: source_term(:, :)
        real (wp),        intent (out)           :: solution(:, :)
        real (wp),        intent (out), optional :: perturbation_optional
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (wp):: perturbation
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Set (real) scalar spherica harmonic coefficients
        !--------------------------------------------------------------------------------

        call this%perform_scalar_analysis( source_term )

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            nlat   => this%NLAT, &
            nlon   => this%NLON, &
            isym   => this%SCALAR_SYMMETRIES, &
            nt     => this%NUMBER_OF_SYNTHESES, &
            xlmbda => helmholtz_constant, &
            sf     => solution, &
            ids    => size( solution, dim = 1), &
            jds    => size( solution, dim = 2), &
            a      => this%workspace%real_harmonic_coefficients, &
            b      => this%workspace%imaginary_harmonic_coefficients, &
            mdab   => size( this%workspace%real_harmonic_coefficients, dim = 1 ), &
            ndab   => size( this%workspace%real_harmonic_coefficients, dim = 2 ), &
            wshsgs => this%workspace%wshsgs, &
            lshsgs => size( this%workspace%wshsgs ), &
            work   => this%workspace%work, &
            lwork  => size( this%workspace%work ), &
            pertrb => perturbation, &
            ierror => error_flag &
            )

            call Islapgs( nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
                mdab, ndab, wshsgs, lshsgs, work, lwork, pertrb, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: invert_HELMHOLTZ'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of VECTOR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for SOLUTION'
                write( stderr, '(A)') 'size( SOLUTION, dim = 1 )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Invalid extent for SOLUTION'
                write( stderr, '(A)') 'size( SOLUTION, dim = 2 )'

            else if ( error_flag == 7 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A) '
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
                write( stderr, '(A)') 'size( SOLUTION, dim = 1 )'

            else if ( error_flag == 8 ) then

                write( stderr, '(A)') 'Invalid extent for '
                write( stderr, '(A)') 'REAL_HARMONIC_COEFFICIENTS (A) '
                write( stderr, '(A)') 'or'
                write( stderr, '(A)') 'IMAGINARY_HARMONIC_COEFFICIENTS (B)'
                write( stderr, '(A)') 'size( SOLUTION, dim = 2 )'

            else if ( error_flag == 9 ) then

                write( stderr, '(A)') 'Invalid extent for WSHSGS'
                write( stderr, '(A)') 'size( WSHSGS )'

            else if ( error_flag == 10 ) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Address optional arguments
        !--------------------------------------------------------------------------------

        if ( present( perturbation_optional ) ) then

            perturbation_optional = perturbation

        end if

    end subroutine invert_helmholtz
    !
    
    ! TODO
    subroutine get_vector_laplacian( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine get_vector_laplacian
    !
    
    ! TODO
    subroutine invert_vector_laplacian( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine invert_vector_laplacian
    !
    
    ! TODO
    subroutine get_stream_function_and_velocity_potential( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine get_stream_function_and_velocity_potential
    !
    
    ! TODO
    subroutine invert_stream_function_and_velocity_potential( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine invert_stream_function_and_velocity_potential
    !
    
    ! TODO
    subroutine perform_grid_transfers( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine perform_grid_transfers
    !
    
    ! TODO
    subroutine perform_geo_math_coordinate_transfers( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine perform_geo_math_coordinate_transfers
    !
    
    !
    subroutine perform_scalar_analysis( this, scalar_function )
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
        class (SpherepackWrapper), intent (in out)   :: this
        real (wp),        intent (in)       :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

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

            call shags( nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                wshags, lshags, work, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag == 0 ) then

            return

        else
            write( stderr, '(A)') 'SPHEREPACK 3.2 error: perform_SCALAR_ANALYSIS'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Invalid extent for WSHAGS'
                write( stderr, '(A)') 'size( WSHAGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for DWORK'
                write( stderr, '(A)') 'size( DWORK )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Failure in GAQD to compute gaussian points '
                write( stderr, '(A)') '(due to failure in eigenvalue routine)'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine perform_scalar_analysis
    !
    
    !
    subroutine perform_scalar_synthesis( this, scalar_function )
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
        class (SpherepackWrapper), intent (in out)  :: this
        real (wp),        intent (out)     :: scalar_function(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
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

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'SPHEREPACK 3.2 error: perform_SCALAR_synthESIS'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Invalid extent for WSHSGS'
                write( stderr, '(A)') 'size( WSHSGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Invalid extent for WORK'
                write( stderr, '(A)') 'size( WORK )'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Invalid extent for DWORK'
                write( stderr, '(A)') 'size( DWORK )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Failure in GAQD to compute gaussian points '
                write( stderr, '(A)') '(due to failure in eigenvalue routine)'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine perform_scalar_synthesis
    !
    
    ! TODO
    subroutine perform_scalar_projection( this, scalar_function, scalar_projection )
        !
        !< Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shpg.html
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out)        :: this
        real (wp), dimension (:, :), intent (in)  :: scalar_function
        real (wp), dimension (:, :), intent (out) :: scalar_projection
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

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

    end subroutine perform_scalar_projection
    !
    
    !
    subroutine perform_vector_analysis( this, vector_function )
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
        class (SpherepackWrapper), intent (in out)  :: this
        real (wp),        intent (in)      :: vector_function(:, :, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        real (wp), allocatable :: polar_component(:, :)
        real (wp), allocatable :: azimuthal_component(:, :)
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            ! Allocate arrays
            allocate( &
                polar_component(     1:nlat, 1:nlon), &
                azimuthal_component( 1:nlat, 1:nlon), &
                stat=allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
                write( stderr, '(A)' ) 'Allocation failed in perform_VECTOR_ANALYSIS'
                write( stderr, '(A)' ) trim( error_message )

            end if

        end associate

        !--------------------------------------------------------------------------------
        ! compute the spherical angle components
        !--------------------------------------------------------------------------------

        associate( &
            F       => vector_function, &
            F_theta => polar_component, &
            F_phi   => azimuthal_component &
            )

            call this%unit_vectors%get_spherical_angle_components( f, f_theta, f_phi )

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

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'ERROR: get_VORTICITY'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of VECTOR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'

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
        ! Release memory
        !--------------------------------------------------------------------------------

        deallocate( &
            polar_component, &
            azimuthal_component, &
            stat=allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( deallocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (SpherepackWrapper)'
            write( stderr, '(A)' ) 'Deallocation failed in PERFORM_VECTOR_ANALYSIS'
            write( stderr, '(A)' ) trim( error_message )

        end if


    end subroutine perform_vector_analysis
    !
    
    !
    subroutine perform_vector_synthesis( this, polar_component, azimuthal_component )
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
        class (SpherepackWrapper), intent (in out)   :: this
        real (wp),                 intent (out)      :: polar_component(:, :)
        real (wp),                 intent (out)      :: azimuthal_component(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: error_flag
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

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

            call vhsgs( nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag == 0 ) then

            return

        else

            write( stderr, '(A)') 'ERROR: get_VORTICITY'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'


            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of VECTOR_SYMMETRIES'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of NUMBER_OF_synthESES'

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

    end subroutine perform_vector_synthesis
    !
    
    ! TODO
    subroutine get_legendre_functions( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine get_legendre_functions
    !
    
    ! TODO
    subroutine Icosahedral_geodesic( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Icosahedral_geodesic
    !
    
    ! TODO
    subroutine perform_multiple_ffts( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%assert_initialized()

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine perform_multiple_ffts
    !
    
    !
    subroutine finalize_spherepack_wrapper( this )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (SpherepackWrapper), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_spherepack_wrapper
    !
    
    !
end module type_SpherepackWrapper
