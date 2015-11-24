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

    ! Local variables confined to the module
    integer, parameter    :: WP = REAL64    !! 64 bit real
    integer, parameter    :: IP = INT32      !! 32 bit integer
    character(200)        :: error_message   !! Probably long enough
    integer (IP)          :: allocate_status
    integer (IP)          :: deallocate_status
    
    ! Declare the sphere derived data type
    type, public ::  sphere_t

        integer (IP)                         :: nlon = 0       !! number of longitudinal points
        integer (IP)                         :: nlat = 0       !! number of latitudinal points
        integer (IP)                         :: ntrunc = 0     !! triangular truncation limit
        integer (IP)                         :: isym = 0       !! symmetries about the equator for scalar calculations
        integer (IP)                         :: itype = 0      !! symmetries about the equator for vector calculations
        complex (WP), dimension (:), pointer :: spec => null() !! Complex (scalar) coefficients
        type (workspace_t)                   :: workspace
        type (grid_t)                        :: grid
        
        ! Commonly used trigonometric functions
        real (WP), dimension (:), allocatable :: sint           !! sin(theta)
        real (WP), dimension (:), allocatable :: cost           !! cos(theta)
        real (WP), dimension (:), allocatable :: sinp           !! sin(phi)
        real (WP), dimension (:), allocatable :: cosp           !! cos(phi)

        ! The spherical unit vectors
        ! Reference: Arfken (1985, p. 102)
        real (WP), dimension (:,:,:), allocatable :: radial_unit_vector
        real (WP), dimension (:,:,:), allocatable :: polar_unit_vector
        real (WP), dimension (:,:,:), allocatable :: azimuthal_unit_vector

        logical, private :: initialized = .false.

    contains
        
        ! SPHEREPACK 3.2 methods
        !        procedure :: Get_colatitude_derivative !! Vtsgs
        procedure :: Get_Gradient !! Gradgs
        !        procedure :: Invert_gradient !!  Igradgs
        procedure :: Get_Divergence !! Divgs
        !        procedure :: Invert_divergence !!Idivgs
        procedure :: Get_Vorticity !! Vrtgs
        !        procedure :: Invert_vorticity !! Ivrtgs
        !        procedure :: Invert_divg_and_vort !! Idvtgs
        !        procedure :: Get_Scalar_laplacian !! Slapgs
        procedure :: Invert_helmholtz !! Islapgs
        !        procedure :: Get_Vector_laplacian !! Vlapgs
        !        procedure :: Invert_vector_laplacian !! Ivlapgs
        !        procedure :: Get_Stream_function_and_velocity_potential
        !        procedure :: Invert_stream_function_and_velocity_potential
        !        procedure :: Perform_grid_transfers
        !        procedure :: Perform_Geo_math_coordinate_transfers
        procedure :: Perform_scalar_analysis
        procedure :: Perform_scalar_synthesis
        !        procedure :: Get_Perform_scalar_projection !! Shpg
        procedure :: Perform_vector_analysis
        procedure :: Perform_vector_synthesis
        !        procedure :: Get_Legendre_functions
        !        procedure :: Get_Icosahedral_geodesic
        !        procedure :: Get_Multiple_ffts
        !        procedure :: Get_Gaussian_weights_and_points !! Contained in type(grid)!! Known as GAQD in SPHEREPACK 3.2
        !        procedure :: Get_Three_dimensional_sphere_graphics
        
        ! Complex methods
        procedure                   :: Perform_complex_analysis
        procedure                   :: Perform_complex_synthesis

        ! Additional public methods
        procedure, non_overridable  :: Create
        procedure, non_overridable  :: Destroy
        procedure, non_overridable  :: Get_index
        procedure, non_overridable  :: Get_coefficient
        procedure                   :: Compute_surface_integral
        procedure                   :: Get_rotation_operator

    end type sphere_t

contains
    !
    !*****************************************************************************************
    !
    subroutine Create( me, nlat, nlon, isym, grid_type )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: me
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
        if ( me%initialized ) then
            print *, 'ERROR: You must destroy "sphere" before re-instantiating'
            return
        end if

        ! Set constants
        me%nlat = nlat

        me%nlon = nlon

        me%ntrunc = nlat - 1 !! Set triangular truncation

        SIZE_SPEC = nlat * (nlat + 1)/2
        
        ! Allocate pointer for complex spectral coefficients
        allocate ( &
            me%spec( 1:size_spec ), &
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

            call Set_scalar_symmetries( me, isym )

        end if

        ! Set grid and workspace
        if ( present( grid_type) ) then

            call me%grid%Create( nlat, nlon, grid_type )
            call me%workspace%Create( nlat, nlon, grid_type )

        else

            call me%grid%Create( nlat, nlon )
            call me%workspace%Create( nlat, nlon )

        end if

        ! Set frequently used trigonometric functions
        call Set_trigonometric_functions( &
            me, &
            me%grid%latitudes, &
            me%grid%longitudes )

        ! Set spherical unit vectors to compute polar and azimuthal components for vector functions
        call Set_spherical_unit_vectors( &
            me, &
            me%sint, &
            me%cost, &
            me%sinp, &
            me%cosp )

        ! Set status
        me%initialized = .true.
        
    end subroutine Create
    !
    !*****************************************************************************************
    !
    subroutine Destroy( me )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. me%initialized ) return

        ! Deallocate pointer
        if ( associated( me%spec ) ) then
            deallocate( &
                me%spec, &
                stat = deallocate_status, &
                errmsg = error_message )
            if ( deallocate_status /= 0 ) then
                print *, 'Pointer deallocation failed in '&
                    &'destruction of sphere_t object: ',&
                    trim( error_message )
                return
            end if

            nullify( me%spec )

        end if

        if ( allocated( me%sint )) deallocate ( me%sint )
        if ( allocated( me%cost )) deallocate ( me%cost )
        if ( allocated( me%sinp )) deallocate ( me%sinp )
        if ( allocated( me%cosp )) deallocate ( me%cosp )
        if ( allocated( me%radial_unit_vector ))  deallocate ( me%radial_unit_vector )
        if ( allocated( me%polar_unit_vector )) deallocate ( me%polar_unit_vector )
        if ( allocated( me%azimuthal_unit_vector )) deallocate ( me%azimuthal_unit_vector )
        
        call me%grid%Destroy()

        call me%workspace%Destroy()

        ! Reset status
        me%initialized = .false.

    end subroutine Destroy
    !
    !*****************************************************************************************
    !
    function Get_index( me, n, m ) result( return_value )
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
        integer (IP)               :: return_value
        class (sphere_t)             :: me
        integer (IP), intent (in)  :: n, m
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: i, ntrunc
        !--------------------------------------------------------------------------------

        ! Set constant
        ntrunc = me%ntrunc

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
    function Get_coefficient( me, n, m ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)    :: me
        integer (IP), intent (in)            :: n, m
        complex (WP)                         :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ntrunc, nm, nm_conjg
        !--------------------------------------------------------------------------------

        ! set truncation limit
        ntrunc = me%ntrunc

        nm = me%Get_index( n, m )

        nm_conjg = me%Get_index(n, -m )

        if ( m < 0 .and. nm_conjg > 0 ) then

            return_value = &
                ( (-1.0_WP)**(-m) ) &
                * conjg( me%spec(nm_conjg) )

        else if ( nm > 0 ) then

            return_value = me%spec(nm)

        else

            return_value = 0.0_WP

        end if

    end function Get_coefficient
    !
    !*****************************************************************************************
    !
    subroutine Set_scalar_symmetries( me, isym )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
        integer (IP), intent (in)         :: isym
        !--------------------------------------------------------------------------------

        if ( isym == 2 ) then

            me%isym = isym

        else if ( isym == 1) then

            me%isym = isym

        else if ( isym == 0 ) then

            me%isym = isym

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
    subroutine Set_vector_symmetries( me, itype )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
        integer (IP), intent (in)         :: itype
        !--------------------------------------------------------------------------------

        if ( itype == 8 ) then

            me%itype = itype

        else if ( itype == 7) then

            me%itype = itype

        else if ( itype == 6) then

            me%itype = itype

        else if ( itype == 5) then

            me%itype = itype

        else if ( itype == 4) then

            me%itype = itype

        else if ( itype == 3) then

            me%itype = itype

        else if ( itype == 2) then

            me%itype = itype

        else if ( itype == 1) then

            me%itype = itype


        else if ( itype == 0 ) then

            me%itype = itype

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
    subroutine Set_trigonometric_functions( me, theta, phi )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)       :: me
        real (WP), dimension (:), intent (in)   :: theta, phi
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)   ::  nlat, nlon
        !--------------------------------------------------------------------------------

        ! Set constants
        nlat = size( theta )
        nlon = size( phi )

        ! Allocate arrays
        allocate ( &
            me%sint( 1:nlat ), &
            me%cost( 1:nlat ), &
            me%sinp( 1:nlon ), &
            me%cosp( 1:nlon ), &
            stat = allocate_status, &
            errmsg = error_message )
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed in '&
                &'Set_trigonometric_functions: ',&
                trim( error_message )
            return
        end if

        ! Compute trigonometric functions
        me%sint = sin( theta )
        me%cost = cos( theta )
        me%sinp = sin( phi )
        me%cosp = cos( phi )

    end subroutine Set_trigonometric_functions
    !
    !*****************************************************************************************
    !
    subroutine Set_spherical_unit_vectors( &
        me, sint, cost, sinp, cosp )
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
        class (sphere_t), intent (in out)      :: me
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
            me%radial_unit_vector( 1:3, 1:nlat, 1:nlon ), &
            me%polar_unit_vector( 1:3, 1:nlat, 1:nlon ), &
            me%azimuthal_unit_vector( 1:3, 1:nlat, 1:nlon ), &
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
                me%radial_unit_vector(:,k,l) = &
                    [ sint(k) * cosp(l), &
                    sint(k) * sinp(l), &
                    cost(k) ]

                ! set polar unit vector
                me%polar_unit_vector(:,k,l) = &
                    [ cost(k) * cosp(l), &
                    cost(k) * sinp(l), &
                    -sint(k) ]

                ! set azimuthal unit vector
                me%azimuthal_unit_vector(:,k,l) = &
                    [ -sinp(l), &
                    cosp(l), &
                    0.0_WP ]
            end do
        end do

    end subroutine Set_spherical_unit_vectors
    !
    !*****************************************************************************************
    !
    subroutine Perform_complex_analysis( me, scalar_function )
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
        class (sphere_t), intent (in out)        :: me
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ntrunc !< Triangular truncation
        integer (IP) :: m,  n  !< counters
        !--------------------------------------------------------------------------------
        
        ntrunc = me%ntrunc
        
        ! compute the (real) spherical harmonic coefficients
        call me%Perform_scalar_analysis( scalar_function )
        
        ! fill complex array dataspec with result
        me%spec = cmplx( &
            0.5_WP * [((me%workspace%a(m + 1,n + 1), n = m,ntrunc), m = 0,ntrunc)], &
            0.5_WP * [((me%workspace%b(m + 1,n + 1), n = m,ntrunc), m = 0,ntrunc)] &
            , WP)
 
    end subroutine Perform_complex_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_complex_synthesis( me, scalar_function )
        ! 
        ! Purpose: 
        ! converts gridded input array (datagrid) to (complex) spherical harmonic coefficients
        ! (dataspec).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)         :: me
        real (WP), dimension (:,:), intent (out)  :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: ntrunc
        integer (IP) :: m, n, nm
        !--------------------------------------------------------------------------------

        ! set the truncation limit
        ntrunc = me%ntrunc
 
        ! Fill real arrays with contents of spec
        do m = 0, ntrunc
            do n = m, ntrunc
                
                ! set the spectral index
                nm = me%Get_index( n, m )
                
                ! set the real component
                me%workspace%a( m + 1,n + 1 ) = &
                    2.0_WP * real( me%spec(nm), WP )
                
                ! set the imaginary component
                me%workspace%b( m + 1,n + 1 ) = &
                    2.0_WP * aimag( me%spec(nm) )
                
            end do
        end do
        
        ! synthesise the scalar function from the (real) harmonic coeffiients
        call me%Perform_scalar_synthesis( scalar_function )
 
    end subroutine Perform_complex_synthesis
    !
    !*****************************************************************************************
    !
    subroutine Synthesize_from_spec( me, spec, scalar_function )
        ! 
        ! Purpose: 
        ! used mainly for testing the spectral method module
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: me
        complex (WP), dimension (:), intent (in)   :: spec
        real (WP), dimension (:,:), intent (out)   :: scalar_function
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ntrunc
        integer (IP):: m
        integer (IP):: n
        integer (IP):: nm
        !--------------------------------------------------------------------------------

        ! set the truncation limit
        ntrunc = me%ntrunc
 
        ! fill real arrays with contents of spec
        do m = 0, ntrunc
            do n = m, ntrunc
                
                ! set the spectral index
                nm = me%Get_index( n, m )
                
                ! set the real component
                me%workspace%a(m + 1,n + 1) = &
                    2.0_WP * real( spec(nm) )
                
                ! set the imaginary component
                me%workspace%b(m + 1,n + 1) = &
                    2.0_WP * aimag( spec(nm) )
                
            end do
        end do
        
        ! synthesise the scalar function from the (real) harmonic coeffiients
        call me%Perform_scalar_synthesis( scalar_function )
 
    end subroutine Synthesize_from_spec
    !
    !*****************************************************************************************
    !
    subroutine Compute_spherical_angle_components( me, &
        vector_function, &
        polar_component, &
        azimuthal_component )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: me
        real (WP), dimension (:,:,:), intent (in)  :: vector_function
        real (WP), dimension (:,:), intent (out)   :: polar_component
        real (WP), dimension (:,:), intent (out)   :: azimuthal_component
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: nlat
        integer (IP)    :: nlon
        integer (IP)    :: k
        integer (IP)    :: l
        type (vector_t) :: theta
        type (vector_t) :: phi
        type (vector_t) :: vector_field
        !--------------------------------------------------------------------------------
        
        ! Set constants
        nlat = me%nlat
        nlon = me%nlon
        
        ! initialize arrays
        polar_component = 0.0_WP
        azimuthal_component = 0.0_WP
        
        ! calculate the spherical angle components
        do l = 1, nlon
            do k = 1, nlat
                
                ! set the vector function
                vector_field = vector_function(:,k,l)

                ! set the latitudinal spherical unit vector
                theta = me%polar_unit_vector(:,k,l)

                ! set the longitudinal spherical unit vector
                phi = me%azimuthal_unit_vector(:,k,l)

                ! set the theta component
                polar_component(k,l) = &
                    theta.dot.vector_field
               
                ! set the azimuthal_component
                azimuthal_component(k,l) = &
                    phi.dot.vector_field

            end do
        end do

    end subroutine Compute_spherical_angle_components
    !
    !*****************************************************************************************
    !

    subroutine Get_rotation_operator( me, &
        scalar_function, rotation_operator)
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: me
        real (WP), dimension (:,:), intent (in)    :: scalar_function
        real (WP), dimension (:,:,:), intent (out) :: rotation_operator
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                            :: nlat
        integer (IP)                            :: nlon
        integer (IP)                            :: k, l !! Counters
        type (vector_t)                         :: theta,  phi
        real (WP), dimension (:,:), allocatable :: polar_gradient_component
        real (WP), dimension (:,:), allocatable :: azimuthal_gradient_component
        !--------------------------------------------------------------------------------

        ! Set constants
        nlat = me%nlat
        nlon = me%nlon

        ! Allocate arrays
        allocate ( &
            polar_gradient_component( 1:nlat, 1:nlon ), &
            azimuthal_gradient_component( 1:nlat, 1:nlon ), &
            stat = allocate_status )
        if ( allocate_status /= 0 ) then
            print *, "Allocation failed!"
            return
        end if

        ! calculate the spherical surface gradient components
        call me%Get_gradient( &
            scalar_function, &
            polar_gradient_component, &
            azimuthal_gradient_component)

        ! initialize array
        rotation_operator = 0.0_WP

        ! calculate the rotation operator applied to a scalar function
        do l = 1, nlon
            do k = 1, nlat

                ! set the theta spherical unit vector
                theta = me%polar_unit_vector(:,k,l)

                ! set the phi spherical unit vector
                phi = me%azimuthal_unit_vector(:,k,l)

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
    function Compute_surface_integral( me, scalar_function ) result( return_value )
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
        class (sphere_t), intent (in out)        :: me
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        real (WP)                                :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                   :: k         !! counter
        real (WP), dimension (me%nlat) :: integrant !! integrant
        !--------------------------------------------------------------------------------

        ! Initialize array
        integrant = 0.0_WP

        ! Compute the integrant
        do k = 1, size(integrant)

            integrant(k) = &
                sum( scalar_function(k,:) ) &
                * me%grid%mesh_phi

        end do

        integrant = &
            me%grid%gaussian_weights * integrant

        return_value = sum( integrant )

    end function Compute_surface_integral
    !
    !*****************************************************************************************
    !
    subroutine Compute_first_moment( me, scalar_function, first_moment )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: me
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        type (vector_t), intent (out)            :: first_moment
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: nlat
        integer (IP)    :: nlon
        integer (IP)    :: k,  l
        type (vector_t) :: u
        real (WP), dimension (me%nlat, me%nlon, 3) :: integrant
        !--------------------------------------------------------------------------------

        ! Set constants
        nlat = me%nlat
        nlon = me%nlon

        ! Initialize array
        integrant = 0.0_WP

        ! Compute integrant
        do l = 1, nlon
            do k = 1, nlat

                u = me%radial_unit_vector(:, k,l)

                integrant(k,l,1) = u%x * scalar_function(k,l)
                integrant(k,l,2) = u%y * scalar_function(k,l)
                integrant(k,l,3) = u%z * scalar_function(k,l)

            end do
        end do

        ! set first moment
        first_moment%x = &
            me%Compute_surface_integral( integrant(:,:, 1))

        first_moment%y = &
            me%Compute_surface_integral( integrant(:,:, 2))

        first_moment%z = &
            me%Compute_surface_integral( integrant(:,:, 3))

    end subroutine Compute_first_moment
    !
    !*****************************************************************************************
    !
    ! SPHEREPACK 3.2 methods
    !
    !*****************************************************************************************
    !
    subroutine Get_colatitude_derivative( me, polar_component, azimuthal_component )
        !
        ! Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vtsgs.html
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: me
        real (WP), dimension (:,:), intent (out) :: polar_component     !! vt
        real (WP), dimension (:,:), intent (out) :: azimuthal_component !! wt
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        !        ! TODO incorporate Vtsgsi into type(workspace)
        !        subroutine vtsgsi(nlat,nlon,wvts,lwvts,work,lwork,dwork,ldwork, ierror)

        !        call Vtsgs( &
        !            size( me%grid%latitudes ), size( me%grid%longitudes ),  me%ityp, 1, &
        !            polar_component, azimuthal_component, &
        !            me%nlat, me%nlon, &
        !            me%workspace%br, me%workspace%bi, &
        !            me%workspace%cr, me%workspace%ci, &
        !            size(me%workspace%br, dim = 1), size(me%workspace%br, dim = 2), &
        !            me%workspace%wvts, me%workspace%size(wvts), &
        !            me%workspace%work, size(me%workspace%work), ierror)

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
    subroutine Get_gradient(  me, scalar_function, &
        polar_gradient_component, &
        azimuthal_gradient_component )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: me
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        real (WP), dimension (:,:), intent (out) :: polar_gradient_component
        real (WP), dimension (:,:), intent (out) :: azimuthal_gradient_component
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------
        
        ! compute the (real) harmonic coefficients
        call me%Perform_scalar_analysis( scalar_function )

        ! compute surface gradient components using coefficients a, b
        !
        ! see: https://www2.cisl.ucar.edu/spherepack/documentation#gradgs.html
        call Gradgs( &
            me%nlat, me%nlon, &
            me%isym, 1, polar_gradient_component, &
            azimuthal_gradient_component, &
            me%nlat, me%nlon, &
            me%workspace%a, me%workspace%b, &
            me%nlat, me%nlat, &
            me%workspace%wvhsgs, &
            size( me%workspace%wvhsgs ), me%workspace%work, &
            size( me%workspace%work ), ierror)

        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Gradgs'
            return
        end if

    end subroutine Get_gradient
    !
    !*****************************************************************************************
    !
    subroutine Invert_gradient( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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

    end subroutine Invert_gradient
    !
    !*****************************************************************************************
    !
    subroutine Get_divergence( me, vector_field, divergence )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: me
        real (WP), dimension (:,:,:), intent (in)  :: vector_field
        real (WP), dimension (:,:), intent (out)   :: divergence
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! calculate the (real) vector harmonic coefficients
        call me%Perform_vector_analysis( vector_field )

        ! calculate the surface divergence
        call Divgs( &
            me%nlat, me%nlon, &
            me%isym, 1, divergence, me%nlat, me%nlon, &
            me%workspace%br, me%workspace%bi, &
            me%nlat, me%nlat, &
            me%workspace%wshsgs, size( me%workspace%wshsgs ), &
            me%workspace%work, size( me%workspace%work ), ierror)

        ! check the error flag
        if (ierror  /=  0) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, 'in Divgs'
            return
        end if

    end subroutine Get_divergence
    !
    !*****************************************************************************************
    !
    subroutine Invert_divergence( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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

    end subroutine Invert_divergence
    !
    !*****************************************************************************************
    !
    subroutine Get_vorticity( me, vector_field, vorticity )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: me
        real (WP), dimension (:,:,:), intent (in)  :: vector_field
        real (WP), dimension (:,:), intent (out)   :: vorticity
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! calculate the (real) vector harmonic coefficients
        call me%Perform_vector_analysis( vector_field )

        ! calculate the surface vorticity
        call Vrtgs( &
            me%nlat,  &
            me%nlon, &
            me%isym, 1, vorticity, &
            me%nlat, me%nlon, &
            me%workspace%cr, me%workspace%ci, &
            me%nlat, me%nlon, &
            me%workspace%wshsgs, size( me%workspace%wshsgs ), &
            me%workspace%work, size( me%workspace%work ), ierror)

        ! check the error flag
        if (ierror  /=  0)  then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Vrtgs'
            return
        end if

    end subroutine Get_vorticity
    !
    !*****************************************************************************************
    !
    subroutine Invert_vorticity( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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

    end subroutine Invert_vorticity
    !
    !*****************************************************************************************
    !
    subroutine Invert_divergence_and_vorticity( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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

    end subroutine Invert_divergence_and_vorticity
    !
    !*****************************************************************************************
    !
    subroutine Get_scalar_laplacian( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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

    end subroutine Get_scalar_laplacian
    !
    !*****************************************************************************************
    !
    subroutine Invert_helmholtz( me, helmholtz_constant, source_term, solution )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)         :: me
        real (WP), intent (in)                    :: helmholtz_constant
        real (WP), dimension (:,:), intent (in)   :: source_term
        real (WP), dimension (:,:), intent (out)  :: solution
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (WP)     :: perturbation
        integer (IP)  :: ierror
        !--------------------------------------------------------------------------------

        call me%Perform_scalar_analysis( source_term )

        ! invert the helmholtz (or poisson) equation
        !
        ! see: https://www2.cisl.ucar.edu/spherepack/documentation#islapgs.html
        call Islapgs( &
            me%nlat, me%nlon, &
            me%isym, 1, helmholtz_constant, &
            solution, me%nlat, me%nlon, &
            me%workspace%a, me%workspace%b, &
            me%nlat, me%nlat, &
            me%workspace%wshsgs, size( me%workspace%wshsgs ), &
            me%workspace%work, size( me%workspace%work ), &
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
    subroutine Get_vector_laplacian( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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

    end subroutine Get_vector_laplacian
    !
    !*****************************************************************************************
    !
    subroutine Invert_vector_laplacian( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
    subroutine Get_stream_function_and_velocity_potential( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
    subroutine Invert_stream_function_and_velocity_potential( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
    subroutine Perform_grid_transfers( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
    subroutine Perform_geo_math_coordinate_transfers( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
    subroutine Perform_scalar_analysis( me, scalar_function )
        !
        ! Purpose:
        ! Converts gridded input array (datagrid) to real spectral coefficients
        ! (dataspec).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)      :: me
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
            me%nlat, me%nlon, &
            me%isym, 1, scalar_function, &
            me%nlat, me%nlon, &
            me%workspace%a, me%workspace%b, &
            me%nlat, me%nlat, &
            me%workspace%wshags, size( me%workspace%wshags ), &
            me%workspace%work, size( me%workspace%work ), ierror )

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shags'
            return
        end if

    end subroutine Perform_scalar_analysis
    !
    !*****************************************************************************************
    !
    subroutine Perform_scalar_synthesis( me, scalar_function )
        !
        ! Purpose:
        ! Converts (real) spherical harmonic coefficients to gridded data array (datagrid).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)       :: me
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
            me%nlat, me%nlon, &
            me%isym, 1, scalar_function, &
            me%nlat, me%nlon, &
            me%workspace%a, me%workspace%b, &
            me%nlat, me%nlat, &
            me%workspace%wshsgs, size( me%workspace%wshsgs ), &
            me%workspace%work, size( me%workspace%work ), ierror)

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror, ' in Shsgs'
            return
        end if

    end subroutine Perform_scalar_synthesis
    !
    !*****************************************************************************************
    !
    subroutine Perform_scalar_projection( me, scalar_function, scalar_projection )
        !
        ! Purpose:
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#shpg.html
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: me
        real (WP), dimension (:,:), intent (in)  :: scalar_function
        real (WP), dimension (:,:), intent (out) :: scalar_projection
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! TODO: Include driver program into type(workspace)
!        call Shpgi( &
!            me%nlat, me%nlon, me%isym, me%ntrunc, &
!            me%workspace%wshp, size( me%workspace%wshp ), &
!            me%workspace%iwshp, size( me%workspace%iwshp ), &
!            me%workspace%work, size( me%workspace%work ), ierror )
!
!        call Shpg( &
!            me%nlat, me%nlon, me%isym, me%ntrunc, &
!            scalar_function, scalar_projection, me%nlat, &
!            me%workspace%wshp, size( me%workspace%wshp ), &
!            me%workspace%iwshp, size( me%workspace%iwshp ), &
!            me%workspace%work, size( me%workspace%work ), ierror )

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
    subroutine Perform_vector_analysis( me, vector_function )
        !
        ! Purpose:
        ! converts gridded input array (datagrid) to real spectral coefficients
        ! (dataspec).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)          :: me
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
        nlat = me%nlat
        nlon = me%nlon

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
        call Compute_spherical_angle_components( &
            me, &
            vector_function, &
            polar_component, &
            azimuthal_component)

        ! calculate (real) vector spherical harmonic analysis
        call Vhags( nlat, nlon, &
            me%isym, 1, polar_component, &
            azimuthal_component, &
            nlat, nlon, me%workspace%br, me%workspace%bi, &
            me%workspace%cr, me%workspace%ci, nlat, nlat, &
            me%workspace%wvhags, size( me%workspace%wvhags ), &
            me%workspace%work, size( me%workspace%work ), ierror)

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
    subroutine Perform_vector_synthesis( me, polar_component, azimuthal_component )
        ! Purpose:
        ! converts gridded input array (datagrid) to real spectral coefficients
        ! (dataspec).
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out)        :: me
        real (WP), dimension (:,:), intent (out) :: polar_component
        real (WP), dimension (:,:), intent (out) :: azimuthal_component
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP):: ierror
        !--------------------------------------------------------------------------------

        ! resynthesise the components from the (real) vector coefficients
        call Vhsgs( &
            me%nlat, me%nlon, &
            me%isym, 1, &
            polar_component,&
            azimuthal_component, &
            me%nlat, me%nlon, &
            me%workspace%br, me%workspace%bi, &
            me%workspace%cr, me%workspace%ci, &
            me%nlat, me%nlat, &
            me%workspace%wvhsgs, size( me%workspace%wvhsgs ), &
            me%workspace%work, size( me%workspace%work ), ierror)

        ! check the error status
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror,' in Vhsgs'
            return
        end if

    end subroutine Perform_vector_synthesis
    !
    !*****************************************************************************************
    !
    subroutine Get_legendre_functions( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
    subroutine Icosahedral_geodesic( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
    subroutine Perform_multiple_ffts( me )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (sphere_t), intent (in out) :: me
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
end module type_sphere_mod
