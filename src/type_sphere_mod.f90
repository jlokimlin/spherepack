!*****************************************************************************************
!
!< Author:
! Jon Lo Kim Lin
!
!< Purpose:
!
! A modern Fortran (2003+) object-oriented spherepack wrapper for NCAR's SPHEREPACK 3.2
!
module type_sphere_mod

    use type_workspace_mod, only: &
        workspace_t

    use type_grid_mod, only: &
        grid_t

    use type_vector_mod
    
    use, intrinsic :: iso_fortran_env, only: &
        REAL64, &
        INT32

    ! Explicit typing only
    implicit none

    ! Everything is private except the derived data type itself
    private
    
    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    integer, parameter    :: WP = REAL64       !! 64 bit real
    integer, parameter    :: IP = INT32        !! 32 bit integer
    character(200)        :: error_message     !! Probably long enough
    integer (IP)          :: allocate_status   !! To check allocation status
    integer (IP)          :: deallocate_status !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, public ::  sphere_t

        ! Components
        logical                               :: initialized = .false. !! Instantiation status
        integer (IP)                          :: nlon = 0              !! number of longitudinal points
        integer (IP)                          :: nlat = 0              !! number of latitudinal points
        integer (IP)                          :: ntrunc = 0            !! triangular truncation limit
        integer (IP)                          :: isym = 0              !! symmetries about the equator for scalar calculations
        integer (IP)                          :: itype = 0             !! symmetries about the equator for vector calculations
        complex (WP), dimension (:), pointer  :: spec => null()        !! Complex (scalar) coefficients
        type (workspace_t), private           :: workspace             !! Contains the various workspace arrays to invoke SPHERPACK 3.2
        type (grid_t)                         :: grid
        
        ! Commonly used trigonometric functions
        real (WP), dimension (:), allocatable :: sint                 !! sin(theta): 0 <= theta <= pi
        real (WP), dimension (:), allocatable :: cost                 !! cos(theta): 0 <= theta <= pi
        real (WP), dimension (:), allocatable :: sinp                 !! sin(phi):   0 <=  phi  <= 2*pi
        real (WP), dimension (:), allocatable :: cosp                 !! cos(phi):   0 <=  phi  <= 2*pi

        ! The spherical unit vectors
        ! Reference: Arfken (1985, p. 102)
        real (WP), dimension (:,:,:), allocatable :: radial_unit_vector
        real (WP), dimension (:,:,:), allocatable :: polar_unit_vector
        real (WP), dimension (:,:,:), allocatable :: azimuthal_unit_vector

    contains
        
        ! All method are private unless stated otherwise
        private

        ! Public SPHEREPACK 3.2 methods
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
        !procedure, public                    :: Get_Gaussian_weights_and_points !! Contained in type(grid)!! Known as GAQD in SPHEREPACK 3.2
        !procedure, public                    :: Get_Three_dimensional_sphere_graphics
        
        ! Complex methods
        procedure, public                    :: Perform_complex_analysis
        procedure, public                    :: Perform_complex_synthesis

        ! Additional public methods
        procedure, non_overridable, public   :: Create
        procedure, non_overridable, public   :: Destroy
        procedure, non_overridable, public   :: Get_index
        procedure, non_overridable, public   :: Get_coefficient
        procedure, public                    :: Compute_surface_integral
        procedure, public                    :: Get_rotation_operator

        ! Private methods
        procedure, non_overridable           :: Assert_initialized
        procedure                            :: Set_trigonometric_functions
        procedure                            :: Set_spherical_unit_vectors
        procedure                            :: Get_spherical_angle_components
        procedure                            :: Set_scalar_symmetries
        procedure                            :: Set_vector_symmetries
        final                                :: Finalize

    end type sphere_t

contains
    !
    !*****************************************************************************************
    !
    subroutine Create( this, nlat, nlon, isym, itype, grid_type )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        integer (IP), intent (in)            :: nlat
        integer (IP), intent (in)            :: nlon
        integer (IP), intent (in), optional  :: isym      !! Either 0, 1, or 2
        integer (IP), intent (in), optional  :: itype     !! Either 0, 1, 2, 3, ..., 8
        character (*), intent (in), optional :: grid_type !! Either '(GAU)' or '(REG)'
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: size_spec
        !--------------------------------------------------------------------------------

        ! Check status
        if ( this%initialized ) then
            print *, 'ERROR: You must destroy "sphere" before re-instantiating'
            return
        end if

        ! Set constants
        this%nlat = nlat

        this%nlon = nlon

        this%ntrunc = nlat - 1 !! Set triangular truncation

        SIZE_SPEC = nlat * (nlat + 1)/2
        
        ! Allocate pointer for complex spectral coefficients
        allocate ( &
            this%spec( 1:size_spec ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            print *, 'Pointer allocation failed in '&
                &'creation of sphere_t object: ',&
                trim( error_message )
            return
        end if

        ! Set scalar symmetries
        if ( present( isym ) ) then

            call this%Set_scalar_symmetries( isym )

        end if

        ! Set vector symmetries
        if (present (itype) ) then

            call this%Set_vector_symmetries( itype )

        end if

        ! Set grid and workspace
        if ( present( grid_type) ) then

            call this%grid%Create( nlat, nlon, grid_type )
            call this%workspace%Create( nlat, nlon, grid_type )

        else

            call this%grid%Create( nlat, nlon )
            call this%workspace%Create( nlat, nlon )

        end if

        ! Set frequently used trigonometric functions
        call this%Set_trigonometric_functions( &
            this%grid%latitudes, &
            this%grid%longitudes )

        ! Set spherical unit vectors to compute polar and azimuthal components for vector functions
        call this%Set_spherical_unit_vectors( &
            this%sint, &
            this%cost, &
            this%sinp, &
            this%cosp )

        ! Set status
        this%initialized = .true.
        
    end subroutine Create
    !
    !*****************************************************************************************
    !
    subroutine Destroy( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. this%initialized ) return

        ! Reset constants
        this%nlon = 0
        this%nlat = 0
        this%ntrunc = 0
        this%isym = 0
        this%itype = 0

        ! Deallocate pointer
        if ( associated( this%spec ) ) then
            deallocate( &
                this%spec, &
                stat = deallocate_status, &
                errmsg = error_message )
            if ( deallocate_status /= 0 ) then
                print *, 'Pointer deallocation failed in '&
                    &'destruction of sphere_t object: ',&
                    trim( error_message )
                return
            end if

            nullify( this%spec )

        end if

        if ( allocated( this%sint )) deallocate ( this%sint )
        if ( allocated( this%cost )) deallocate ( this%cost )
        if ( allocated( this%sinp )) deallocate ( this%sinp )
        if ( allocated( this%cosp )) deallocate ( this%cosp )
        if ( allocated( this%radial_unit_vector ))  deallocate ( this%radial_unit_vector )
        if ( allocated( this%polar_unit_vector )) deallocate ( this%polar_unit_vector )
        if ( allocated( this%azimuthal_unit_vector )) deallocate ( this%azimuthal_unit_vector )
        
        call this%grid%Destroy()

        call this%workspace%Destroy()

        ! Reset status
        this%initialized = .false.

    end subroutine Destroy
    !
    !*****************************************************************************************
    !
    subroutine Assert_initialized( this )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. this%initialized ) then
            print *, 'ERROR: You must instantiate "sphere" object '&
                &//'before calling methods'
            stop
        end if

    end subroutine Assert_initialized
    !
    !*****************************************************************************************
    !
    function Get_index( this, n, m ) result( return_value )
        !
        ! Purpose:
        ! The spectral data is assumed to be in a complex array of dimension
        ! (MTRUNC+1)*(MTRUNC+2)/2. MTRUNC is the triangular truncation limit
        ! (MTRUNC = 42 for T42). MTRUNC must be <= nlat-1.
        !
        ! The coefficients are ordered so that
        !
        ! first (nm=1)  is m=0,n=0, second (nm=2) is m=0,n=1,
        ! nm=MTRUNC is m=0,n=MTRUNC, nm=MTRUNC+1 is m=1,n=1, etc.
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
        ! integer (IP), dimension ((MTRUNC+1)*(MTRUNC+2)/2) :: indxm,indxn
        ! indxm = [((m,n=m,MTRUNC),m=0,MTRUNC)]
        ! indxn = [((n,n=m,MTRUNC),m=0,MTRUNC)]
        !
        ! Conversely, the index nm as a function of m and n is:
        ! nm = sum([(i,i=MTRUNC+1,MTRUNC-m+2,-1)])+n-m+1
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)                       :: return_value
        class (sphere_t), intent (in out)  :: this
        integer (IP), intent (in)          :: n, m
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: i, ntrunc
        !--------------------------------------------------------------------------------

        ! Set constant
        ntrunc = this%ntrunc

        if ( m <= n .and. max(n, m) <= ntrunc ) then

            return_value = &
                sum ( [(i, i = ntrunc+1, ntrunc-m+2, -1)] ) + n-m+1
        else

            return_value = -1

        end if


    end function Get_index
    !
    !*****************************************************************************************
    !
    function Get_coefficient( this, n, m ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: this
        integer (IP), intent (in)            :: n, m
        complex (WP)                         :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ntrunc, nm, nm_conjg
        !--------------------------------------------------------------------------------

        ! set truncation limit
        ntrunc = this%ntrunc

        nm = this%Get_index( n, m )

        nm_conjg = this%Get_index(n, -m )

        if ( m < 0 .and. nm_conjg > 0 ) then

            return_value = &
                ( (-1.0_WP)**(-m) ) &
                * conjg( this%spec(nm_conjg) )

        else if ( nm > 0 ) then

            return_value = this%spec(nm)

        else

            return_value = 0.0_WP

        end if

    end function Get_coefficient
    !
    !*****************************************************************************************
    !
    subroutine Set_scalar_symmetries( this, isym )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        integer (IP), intent (in)         :: isym
        !--------------------------------------------------------------------------------

        if ( isym == 2 ) then

            this%isym = isym

        else if ( isym == 1) then

            this%isym = isym

        else if ( isym == 0 ) then

            this%isym = isym

        else
            ! Handle invalid isym
            print *, 'ERROR: optional argument isym = ', isym
            print *, 'must be either 0, 1, or 2 (default isym = 0)'
            stop
        end if

    end subroutine Set_scalar_symmetries
        !
    !*****************************************************************************************
    !
    subroutine Set_vector_symmetries( this, itype )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: this
        integer (IP), intent (in)         :: itype
        !--------------------------------------------------------------------------------

        if ( itype == 8 ) then

            this%itype = itype

        else if ( itype == 7) then

            this%itype = itype

        else if ( itype == 6) then

            this%itype = itype

        else if ( itype == 5) then

            this%itype = itype

        else if ( itype == 4) then

            this%itype = itype

        else if ( itype == 3) then

            this%itype = itype

        else if ( itype == 2) then

            this%itype = itype

        else if ( itype == 1) then

            this%itype = itype


        else if ( itype == 0 ) then

            this%itype = itype

        else
            ! Handle invalid isym
            print *, 'ERROR: optional argument itype = ', itype
            print *, 'must be either 0, 1, 2, ..., 8 (default itype = 0)'
            stop

        end if

    end subroutine Set_vector_symmetries
    !
    !*****************************************************************************************
    !
    subroutine Set_trigonometric_functions( this, theta, phi )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)       :: this
        real (WP), dimension (:), intent (in)   :: theta, phi
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)   ::  nlat, nlon
        !--------------------------------------------------------------------------------

        ! Check if magnetization is allocated
        if ( .not. allocated( this%grid%latitudes ) ) then
            print *, 'ERROR: You must allocate "latitudes" '&
                &//'before calling "Set_trigonometric_functions"'
            stop
        else if ( .not. allocated( this%grid%longitudes ) ) then
            print *, 'ERROR: You must allocate "longitudes" '&
                &//'before calling "Set_trigonometric_functions"'
            stop
        end if

        ! Set constants
        nlat = size( theta )
        nlon = size( phi )

        ! Allocate arrays
        allocate ( &
            this%sint( 1:nlat ), &
            this%cost( 1:nlat ), &
            this%sinp( 1:nlon ), &
            this%cosp( 1:nlon ), &
            stat = allocate_status, &
            errmsg = error_message )
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &'Set_trigonometric_functions: ',&
                trim( error_message )
            return
        end if

        ! Compute trigonometric functions
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
        ! Purpose:
        ! Sets the spherical unit vectors
        !
        ! Remark:
        ! The "grid" component of sphere must be
        ! initialized before calling this subroutine
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)      :: this
        real (WP), dimension (:), intent (in)  :: sint, cost, sinp, cosp
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  nlat, nlon, k, l
        !--------------------------------------------------------------------------------

        ! Set constants
        nlat = size( sint )
        nlon = size( sinp )

        ! Allocate arrays
        allocate ( &
            this%radial_unit_vector( 1:3, 1:nlat, 1:nlon ), &
            this%polar_unit_vector( 1:3, 1:nlat, 1:nlon ), &
            this%azimuthal_unit_vector( 1:3, 1:nlat, 1:nlon ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &'Set_spherical_unit_vectors: ',&
                trim( error_message )
            return
        end if

        ! Compute spherical unit vectors
        do l = 1, nlon
            do k = 1, nlat

                ! set radial unit vector
                this%radial_unit_vector(:,k,l) = &
                    [ sint(k) * cosp(l), &
                    sint(k) * sinp(l), &
                    cost(k) ]

                ! set polar unit vector
                this%polar_unit_vector(:,k,l) = &
                    [ cost(k) * cosp(l), &
                    cost(k) * sinp(l), &
                    -sint(k) ]

                ! set azimuthal unit vector
                this%azimuthal_unit_vector(:,k,l) = &
                    [ -sinp(l), &
                    cosp(l), &
                    0.0_WP ]
            end do
        end do

    end subroutine Set_spherical_unit_vectors
    !
    !*****************************************************************************************
    !
    subroutine Perform_complex_analysis( this, scalar_function )
        ! 
        ! Purpose: 
        ! converts gridded input array (scalar_function) to (complex) spherical harmonic coefficients
        ! (dataspec).
        !
        ! the spectral data is assumed to be in a complex array of dimension
        ! (mtrunc+1)*(mtrunc+2)/2. mtrunc is the triangular truncation limit
        ! (mtrunc = 42 for t42). mtrunc must be <= nlat-1. coefficients are
        ! ordered so that first (nm=1) is m=0,n=0, second is m=0,n=1, 
        ! nm=mtrunc is m=0,n=mtrunc, nm=mtrunc+1 is m=1,n=1, etc.
        ! in fortran95 syntax, values of m (degree) and n (order) as a function
        ! of the index nm are: 

        ! integer (IP), dimension ((mtrunc+1)*(mtrunc+2)/2) :: indxm,indxn
        ! indxm = [((m,n=m,mtrunc),m=0,mtrunc)]
        ! indxn = [((n,n=m,mtrunc),m=0,mtrunc)]

        ! conversely, the index nm as a function of m and n is: 
        ! nm = sum([(i,i=mtrunc+1,mtrunc-m+2,-1)])+n-m+1
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ntrunc !< Triangular truncation
        integer (IP) :: m,  n  !< counters
        !--------------------------------------------------------------------------------
        
        ! Check status
        call this%Assert_initialized()

        ntrunc = this%ntrunc
        
        ! compute the (real) spherical harmonic coefficients
        call this%Perform_scalar_analysis( scalar_function )
        
        ! fill complex array dataspec with result
        this%spec = cmplx( &
            0.5_WP * [((this%workspace%a(m + 1,n + 1), n = m,ntrunc), m = 0,ntrunc)], &
            0.5_WP * [((this%workspace%b(m + 1,n + 1), n = m,ntrunc), m = 0,ntrunc)] &
            , WP)
 
    end subroutine Perform_complex_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_complex_synthesis( this, scalar_function )
        ! 
        ! Purpose: 
        ! converts gridded input array (datagrid) to (complex) spherical harmonic coefficients
        ! (dataspec).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)         :: this
        real (WP), dimension (:,:), intent (out)  :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ntrunc
        integer (IP) :: m, n, nm
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! set the truncation limit
        ntrunc = this%ntrunc
 
        ! Fill real arrays with contents of spec
        do m = 0, ntrunc
            do n = m, ntrunc
                
                ! set the spectral index
                nm = this%Get_index( n, m )
                
                ! set the real component
                this%workspace%a( m + 1,n + 1 ) = &
                    2.0_WP * real( this%spec(nm), WP )
                
                ! set the imaginary component
                this%workspace%b( m + 1,n + 1 ) = &
                    2.0_WP * aimag( this%spec(nm) )
                
            end do
        end do
        
        ! synthesise the scalar function from the (real) harmonic coeffiients
        call this%Perform_scalar_synthesis( scalar_function )
 
    end subroutine Perform_complex_synthesis
    !
    !*****************************************************************************************
    !
    subroutine Synthesize_from_spec( this, spec, scalar_function )
        ! 
        ! Purpose: 
        ! used mainly for testing the spectral method module
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: this
        complex (WP), dimension (:), intent (in)   :: spec
        real (WP), dimension (:,:), intent (out)   :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ntrunc   !! Triangular truncation limit
        integer (IP):: m, n     !! Counters
        integer (IP):: nm       !! Index
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! set the truncation limit
        ntrunc = this%ntrunc
 
        ! fill real arrays with contents of spec
        do m = 0, ntrunc
            do n = m, ntrunc
                
                ! set the spectral index
                nm = this%Get_index( n, m )
                
                ! set the real component
                this%workspace%a(m + 1,n + 1) = &
                    2.0_WP * real( spec(nm) )
                
                ! set the imaginary component
                this%workspace%b(m + 1,n + 1) = &
                    2.0_WP * aimag( spec(nm) )
                
            end do
        end do
        
        ! synthesise the scalar function from the (real) harmonic coeffiients
        call this%Perform_scalar_synthesis( scalar_function )
 
    end subroutine Synthesize_from_spec
    !
    !*****************************************************************************************
    !
    subroutine Get_spherical_angle_components( this, &
        vector_function, polar_component, azimuthal_component )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: this
        real (WP), dimension (:,:,:), intent (in)  :: vector_function
        real (WP), dimension (:,:), intent (out)   :: polar_component
        real (WP), dimension (:,:), intent (out)   :: azimuthal_component
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: nlat, nlon
        integer (IP)    :: k,  l        !! Counters
        type (vector_t) :: theta, phi   !! Spherical unit vectors
        type (vector_t) :: vector_field
        !--------------------------------------------------------------------------------
        
        ! Check status
        call this%Assert_initialized()

        ! Set constants
        nlat = this%nlat
        nlon = this%nlon
        
        ! initialize arrays
        polar_component = 0.0_WP
        azimuthal_component = 0.0_WP
        
        ! calculate the spherical angle components
        do l = 1, nlon
            do k = 1, nlat
                
                ! set the vector function
                vector_field = vector_function(:,k,l)

                ! set the latitudinal spherical unit vector
                theta = this%polar_unit_vector(:,k,l)

                ! set the longitudinal spherical unit vector
                phi = this%azimuthal_unit_vector(:,k,l)

                ! set the theta component
                polar_component(k,l) = &
                    theta.dot.vector_field
               
                ! set the azimuthal_component
                azimuthal_component(k,l) = &
                    phi.dot.vector_field

            end do
        end do

    end subroutine Get_spherical_angle_components
    !
    !*****************************************************************************************
    !
    subroutine Get_rotation_operator( this, scalar_function, rotation_operator)
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: this
        real (WP), dimension (:,:), intent (in)    :: scalar_function
        real (WP), dimension (:,:,:), intent (out) :: rotation_operator
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                            :: nlat, nlon
        integer (IP)                            :: k, l         !! Counters
        type (vector_t)                         :: theta,  phi
        real (WP), dimension (:,:), allocatable :: polar_gradient_component
        real (WP), dimension (:,:), allocatable :: azimuthal_gradient_component
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! Set constants
        nlat = this%nlat
        nlon = this%nlon

        ! Allocate arrays
        allocate ( &
            polar_gradient_component( 1:nlat, 1:nlon ), &
            azimuthal_gradient_component( 1:nlat, 1:nlon ), &
            stat = allocate_status )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, "Allocation failed!"
            return
        end if

        ! calculate the spherical surface gradient components
        call this%Get_gradient( &
            scalar_function, polar_gradient_component, azimuthal_gradient_component)

        ! initialize array
        rotation_operator = 0.0_WP

        ! calculate the rotation operator applied to a scalar function
        do l = 1, nlon
            do k = 1, nlat

                ! set the theta spherical unit vector
                theta = this%polar_unit_vector(:,k,l)

                ! set the phi spherical unit vector
                phi = this%azimuthal_unit_vector(:,k,l)

                rotation_operator(:, k,l) = &
                    phi * polar_gradient_component(k,l) &
                    - theta * azimuthal_gradient_component(k,l)
            end do

        end do

        ! Deallocate arrays
        deallocate ( &
            polar_gradient_component, &
            azimuthal_gradient_component, &
            stat = deallocate_status )
        if ( deallocate_status /= 0 ) then
            print *, 'Deallocation failed in '&
                &//'Get_rotation_operator'
            return
        end if

    end subroutine Get_rotation_operator
    !
    !*****************************************************************************************
    !
    function Compute_surface_integral( this, scalar_function ) result( return_value )
        !
        ! Purpose:
        ! computes the surface integral on the sphere (trapezoidal rule in phi
        ! and gaussian quadrature in theta)
        !
        !   \int_{s^1} sf(theta, phi) ds
        !
        !    for ds : = ds(theta,phi) = sin(theta) dtheta dphi
        !
        !    for 0 <= theta <= pi and 0 <= phi <= 2*pi
        !
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        real (WP)                                :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                     :: k         !! counter
        real (WP), dimension (this%nlat) :: integrant !! integrant
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! Initialize array
        integrant = 0.0_WP

        ! Compute the integrant
        do k = 1, size(integrant)

            integrant(k) = &
                sum( scalar_function(k,:) ) &
                * this%grid%mesh_phi

        end do

        integrant = &
            this%grid%gaussian_weights * integrant

        return_value = sum( integrant )

    end function Compute_surface_integral
    !
    !*****************************************************************************************
    !
    subroutine Compute_first_moment( this, scalar_function, first_moment )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        type (vector_t), intent (out)            :: first_moment
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: nlat
        integer (IP)    :: nlon
        integer (IP)    :: k,  l
        type (vector_t) :: u
        real (WP), dimension (this%nlat, this%nlon, 3) :: integrant
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! Set constants
        nlat = this%nlat
        nlon = this%nlon

        ! Initialize array
        integrant = 0.0_WP

        ! Compute integrant
        do l = 1, nlon
            do k = 1, nlat

                u = this%radial_unit_vector(:, k,l)

                integrant(k,l,1) = u%x * scalar_function(k,l)
                integrant(k,l,2) = u%y * scalar_function(k,l)
                integrant(k,l,3) = u%z * scalar_function(k,l)

            end do
        end do

        ! set first moment
        first_moment%x = &
            this%Compute_surface_integral( integrant(:,:, 1))

        first_moment%y = &
            this%Compute_surface_integral( integrant(:,:, 2))

        first_moment%z = &
            this%Compute_surface_integral( integrant(:,:, 3))

    end subroutine Compute_first_moment
    !
    !*****************************************************************************************
    !
    ! SPHEREPACK 3.2 methods
    !
    !*****************************************************************************************
    !
    subroutine Get_colatitude_derivative( this, polar_component, azimuthal_component )
        !
        ! Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vtsgs.html
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:,:), intent (out) :: polar_component     !! vt
        real (WP), dimension (:,:), intent (out) :: azimuthal_component !! wt
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        !        ! TODO incorporate Vtsgsi into type(workspace)
        !        subroutine vtsgsi(nlat,nlon,wvts,lwvts,work,lwork,dwork,ldwork, ierror)

        !        call Vtsgs( &
        !            size( this%grid%latitudes ), size( this%grid%longitudes ),  this%ityp, 1, &
        !            polar_component, azimuthal_component, &
        !            this%nlat, this%nlon, &
        !            this%workspace%br, this%workspace%bi, &
        !            this%workspace%cr, this%workspace%ci, &
        !            size(this%workspace%br, dim = 1), size(this%workspace%br, dim = 2), &
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
        polar_gradient_component, &
        azimuthal_gradient_component )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        real (WP), dimension (:,:), intent (out) :: polar_gradient_component
        real (WP), dimension (:,:), intent (out) :: azimuthal_gradient_component
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------
        
        ! Check status
        call this%Assert_initialized()

        ! compute the (real) harmonic coefficients
        call this%Perform_scalar_analysis( scalar_function )

        ! compute surface gradient components using coefficients a, b
        !
        ! see: https://www2.cisl.ucar.edu/spherepack/documentation#gradgs.html
        call Gradgs( &
            this%nlat, this%nlon, &
            this%isym, 1, polar_gradient_component, &
            azimuthal_gradient_component, &
            this%nlat, this%nlon, &
            this%workspace%a, this%workspace%b, &
            this%nlat, this%nlat, &
            this%workspace%wvhsgs, &
            size( this%workspace%wvhsgs ), this%workspace%work, &
            size( this%workspace%work ), ierror)

        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Gradgs'
            return
        end if

    end subroutine Get_gradient
    !
    !*****************************************************************************************
    !
    subroutine Invert_gradient( this )
        !
        ! Purpose:
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

    end subroutine Invert_gradient
    !
    !*****************************************************************************************
    !
    subroutine Get_divergence( this, vector_field, divergence )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: this
        real (WP), dimension (:,:,:), intent (in)  :: vector_field
        real (WP), dimension (:,:), intent (out)   :: divergence
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! calculate the (real) vector harmonic coefficients
        call this%Perform_vector_analysis( vector_field )

        ! calculate the surface divergence
        call Divgs( &
            this%nlat, this%nlon, &
            this%isym, 1, divergence, this%nlat, this%nlon, &
            this%workspace%br, this%workspace%bi, &
            this%nlat, this%nlat, &
            this%workspace%wshsgs, size( this%workspace%wshsgs ), &
            this%workspace%work, size( this%workspace%work ), ierror)

        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, 'in Divgs'
            return
        end if

    end subroutine Get_divergence
    !
    !*****************************************************************************************
    !
    subroutine Invert_divergence( this )
        !
        ! Purpose:
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
    subroutine Get_vorticity( this, vector_field, vorticity )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: this
        real (WP), dimension (:,:,:), intent (in)  :: vector_field
        real (WP), dimension (:,:), intent (out)   :: vorticity
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! Check status
        call this%Assert_initialized()

        ! calculate the (real) vector harmonic coefficients
        call this%Perform_vector_analysis( vector_field )

        ! calculate the surface vorticity
        call Vrtgs( &
            this%nlat,  &
            this%nlon, &
            this%isym, 1, vorticity, &
            this%nlat, this%nlon, &
            this%workspace%cr, this%workspace%ci, &
            this%nlat, this%nlon, &
            this%workspace%wshsgs, size( this%workspace%wshsgs ), &
            this%workspace%work, size( this%workspace%work ), ierror)

        ! check the error flag
        if (ierror  /=  0)  then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vrtgs'
            return
        end if

    end subroutine Get_vorticity
    !
    !*****************************************************************************************
    !
    subroutine Invert_vorticity( this )
        !
        ! Purpose:
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
        ! Purpose:
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
        ! Purpose:
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
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)         :: this
        real (WP), intent (in)                    :: helmholtz_constant
        real (WP), dimension (:,:), intent (in)   :: source_term
        real (WP), dimension (:,:), intent (out)  :: solution
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
            this%nlat, this%nlon, &
            this%isym, 1, helmholtz_constant, &
            solution, this%nlat, this%nlon, &
            this%workspace%a, this%workspace%b, &
            this%nlat, this%nlat, &
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
        ! Purpose:
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
        ! Purpose:
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
        ! Purpose:
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
        ! Purpose:
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
        ! Purpose:
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
        ! Purpose:
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
        ! Purpose:
        ! Converts gridded input array (datagrid) to real spectral coefficients
        ! (dataspec).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)      :: this
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! perform the (real) spherical harmonic analysis
        !
        ! see: https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        call Shags( &
            this%nlat, this%nlon, &
            this%isym, 1, scalar_function, &
            this%nlat, this%nlon, &
            this%workspace%a, this%workspace%b, &
            this%nlat, this%nlat, &
            this%workspace%wshags, size( this%workspace%wshags ), &
            this%workspace%work, size( this%workspace%work ), ierror )

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shags'
            return
        end if

    end subroutine Perform_scalar_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_scalar_synthesis( this, scalar_function )
        !
        ! Purpose:
        ! Converts (real) spherical harmonic coefficients to gridded data array (datagrid).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)       :: this
        real (WP), dimension (:,:), intent (out)  :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! perform (real) spherical harmonic synthesis
        !
        ! see: https://www2.cisl.ucar.edu/spherepack/documentation#shsgs.html
        call Shsgs( &
            this%nlat, this%nlon, &
            this%isym, 1, scalar_function, &
            this%nlat, this%nlon, &
            this%workspace%a, this%workspace%b, &
            this%nlat, this%nlat, &
            this%workspace%wshsgs, size( this%workspace%wshsgs ), &
            this%workspace%work, size( this%workspace%work ), ierror)

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shsgs'
            return
        end if

    end subroutine Perform_scalar_synthesis
    !
    !*****************************************************************************************
    !
    subroutine Perform_scalar_projection( this, scalar_function, scalar_projection )
        !
        ! Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shpg.html
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        real (WP), dimension (:,:), intent (out) :: scalar_projection
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! TODO: Include driver program into type(workspace)
        !        call Shpgi( &
        !            this%nlat, this%nlon, this%isym, this%ntrunc, &
        !            this%workspace%wshp, size( this%workspace%wshp ), &
        !            this%workspace%iwshp, size( this%workspace%iwshp ), &
        !            this%workspace%work, size( this%workspace%work ), ierror )
        !
        !        call Shpg( &
        !            this%nlat, this%nlon, this%isym, this%ntrunc, &
        !            scalar_function, scalar_projection, this%nlat, &
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
        ! Purpose:
        ! converts gridded input array (datagrid) to real spectral coefficients
        ! (dataspec).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: this
        real (WP), dimension (:,:,:), intent (in)  :: vector_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                            :: nlat
        integer (IP)                            :: nlon
        integer (IP)                            :: ierror
        real (WP), dimension (:,:), allocatable :: polar_component
        real (WP), dimension (:,:), allocatable :: azimuthal_component
        !--------------------------------------------------------------------------------

        ! Set contants
        nlat = this%nlat
        nlon = this%nlon

        ! Allocate arrays
        allocate ( &
            polar_component( 1:nlat, 1:nlon), &
            azimuthal_component( 1:nlat, 1:nlon), &
            stat = allocate_status )
        if ( allocate_status /= 0 ) then
            print *, "Allocation failed!"
            return
        end if

        ! compute the spherical angle components
        call this%Get_spherical_angle_components( &
            vector_function, &
            polar_component, &
            azimuthal_component)

        ! calculate (real) vector spherical harmonic analysis
        call Vhags( nlat, nlon, &
            this%isym, 1, polar_component, &
            azimuthal_component, &
            nlat, nlon, this%workspace%br, this%workspace%bi, &
            this%workspace%cr, this%workspace%ci, nlat, nlat, &
            this%workspace%wvhags, size( this%workspace%wvhags ), &
            this%workspace%work, size( this%workspace%work ), ierror)

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror,' in Vhags'
            return
        end if

        ! Deallocate arrays
        deallocate ( &
            polar_component, &
            azimuthal_component, &
            stat = deallocate_status )
        if ( deallocate_status /= 0 ) then
            print *, 'Deallocation failed in '&
                &//'Perform_vector_analysis'
            return
        end if

    end subroutine Perform_vector_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_vector_synthesis( this, polar_component, azimuthal_component )
        ! Purpose:
        ! converts gridded input array (datagrid) to real spectral coefficients
        ! (dataspec).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: this
        real (WP), dimension (:,:), intent (out) :: polar_component
        real (WP), dimension (:,:), intent (out) :: azimuthal_component
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! resynthesise the components from the (real) vector coefficients
        call Vhsgs( &
            this%nlat, this%nlon, &
            this%isym, 1, &
            polar_component,&
            azimuthal_component, &
            this%nlat, this%nlon, &
            this%workspace%br, this%workspace%bi, &
            this%workspace%cr, this%workspace%ci, &
            this%nlat, this%nlat, &
            this%workspace%wvhsgs, size( this%workspace%wvhsgs ), &
            this%workspace%work, size( this%workspace%work ), ierror)

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror,' in Vhsgs'
            return
        end if

    end subroutine Perform_vector_synthesis
    !
    !*****************************************************************************************
    !
    subroutine Get_legendre_functions( this )
        !
        ! Purpose:
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

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Get_legendre_functions
    !
    !*****************************************************************************************
    !
    subroutine Icosahedral_geodesic( this )
        !
        ! Purpose:
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

        ierror = 0
        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vtsgs'
            return
        end if

    end subroutine Icosahedral_geodesic
    !
    !*****************************************************************************************
    !
    subroutine Perform_multiple_ffts( this )
        !
        ! Purpose:
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
    subroutine Finalize( this )
        !
        ! Purpose:
        !< Finalize object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (sphere_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%Destroy()

    end subroutine Finalize
    !
    !*****************************************************************************************
    !
end module type_sphere_mod
