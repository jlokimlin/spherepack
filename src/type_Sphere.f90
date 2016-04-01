module type_Sphere

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_Workspace, only: &
        Workspace

    use type_SphericalGrid, only: &
        SphericalGrid

    use type_SphericalUnitVectors, only: &
        SphericalUnitVectors

    use type_TrigonometricFunctions, only: &
        TrigonometricFunctions

    use type_ThreeDimensionalVector, only: &
        Vector => ThreeDimensionalVector, &
        assignment(=), &
        operator(*)
    
    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: Sphere

    ! Declare derived data type
    type, abstract, public :: Sphere
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                            public  :: initialized = .false. !! Instantiation status
        integer (ip),                       public  :: NUMBER_OF_LONGITUDES = 0   !! number of longitudinal points
        integer (ip),                       public  :: NUMBER_OF_LATITUDES = 0   !! number of latitudinal points
        integer (ip),                       public  :: TRIANGULAR_TRUNCATION_LIMIT = 0 !! triangular truncation limit
        integer (ip),                       public  :: SCALAR_SYMMETRIES = 0 !! symmetries about the equator for scalar calculations
        integer (ip),                       public  :: VECTOR_SYMMETRIES = 0 !! symmetries about the equator for vector calculations
        integer (ip),                       public  :: NUMBER_OF_SYNTHESES = 0
        integer (ip),          allocatable, public  :: INDEX_ORDER_M(:)
        integer (ip),          allocatable, public  :: INDEX_DEGREE_N(:)
        real (wp),                          public  :: RADIUS_OF_SPHERE = 0.0_wp
        real (wp),             allocatable, public  :: vorticity_and_divergence_coefficients(:)
        real (wp),             allocatable, public  :: laplacian_coefficients(:)
        real (wp),             allocatable, public  :: inverse_laplacian_coefficients(:)
        complex (wp),          allocatable, public  :: complex_spectral_coefficients(:)
        class (Workspace),     allocatable, public  :: workspace
        class (SphericalGrid), allocatable, public  :: grid
        type (TrigonometricFunctions),      public  :: trigonometric_functions
        type (SphericalUnitVectors),        public  :: unit_vectors
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure,                              public  :: create_sphere
        procedure,                              public  :: destroy => destroy_sphere
        procedure,                              public  :: destroy_sphere
        procedure,                              public  :: get_index
        procedure,                              public  :: get_coefficient
        procedure,                              public  :: get_scalar_laplacian
        procedure,                              public  :: invert_scalar_laplacian
        generic, public :: get_laplacian => get_scalar_laplacian
        generic, public :: invert_laplacian => invert_scalar_laplacian
        procedure,                              public  :: invert_helmholtz
        procedure,                              public  :: get_gradient
        procedure,                              public  :: invert_gradient
        procedure,                              public  :: get_vorticity
        procedure,                              public  :: invert_vorticity
        procedure, private :: get_divergence_from_vector_field
        procedure, private :: get_divergence_from_spherical_angles
        generic,   public   :: get_divergence => &
            get_divergence_from_vector_field, &
            get_divergence_from_spherical_angles
        procedure,                              public  :: invert_divergence
        procedure,                              public  :: get_rotation_operator => compute_angular_momentum
        procedure,                              private :: get_scalar_symmetries
        procedure,                              private :: get_vector_symmetries
        procedure,                              public  :: perform_complex_analysis
        procedure,                              public  :: perform_complex_synthesis
        procedure (scalar_analysis),  deferred, public  :: perform_scalar_analysis
        procedure (scalar_synthesis), deferred, public  :: perform_scalar_synthesis
        procedure,                              public  :: perform_vector_analysis
        procedure (vector_analysis),  deferred, public  :: vector_analysis_from_spherical_components
        procedure (vector_synthesis), deferred, public  :: perform_vector_synthesis
        procedure,                              public  :: synthesize_from_complex_spectral_coefficients
        procedure,                              public  :: analyze_into_complex_spectral_coefficients
        procedure,                              public  :: get_vorticity_and_divergence_coefficients_from_velocities
        procedure,                              public  :: get_velocities_from_vorticity_and_divergence_coefficients
        !----------------------------------------------------------------------
    end type Sphere

    abstract interface
        subroutine scalar_analysis(this, scalar_function)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            class (Sphere), intent (in out) :: this
            real (wp),      intent (in)     :: scalar_function(:, :)
            !----------------------------------------------------------------------
        end subroutine scalar_analysis
    end interface


    abstract interface
        subroutine scalar_synthesis(this, scalar_function)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            class (Sphere), intent (in out) :: this
            real (wp),      intent (out)    :: scalar_function(:, :)
            !----------------------------------------------------------------------
        end subroutine scalar_synthesis
    end interface

    abstract interface
        subroutine vector_analysis(this, &
            polar_component, azimuthal_component)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            class (Sphere), intent (in out) :: this
            real (wp),      intent (in)     :: polar_component(:, :)
            real (wp),      intent (in)     :: azimuthal_component(:, :)
            !----------------------------------------------------------------------
        end subroutine vector_analysis
    end interface

    abstract interface
        subroutine vector_synthesis(this, polar_component, azimuthal_component)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            class (Sphere), intent (in out) :: this
            real (wp),      intent (out)    :: polar_component(:, :)
            real (wp),      intent (out)    :: azimuthal_component(:, :)
            !----------------------------------------------------------------------
        end subroutine vector_synthesis
    end interface

!    ! Overload class method
!    interface get_divergence_interface
!        module procedure get_divergence_from_vector_field
!        module procedure get_divergence_from_spherical_angles
!    end interface get_divergence_interface

contains


    subroutine create_sphere( this, nlat, nlon, ntrunc, isym, itype, isynt, rsphere )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        integer (ip),   intent (in)     :: nlat
        integer (ip),   intent (in)     :: nlon
        integer (ip),   intent (in)     :: ntrunc
        integer (ip),   intent (in)     :: isym  !! Either 0, 1, or 2
        integer (ip),   intent (in)     :: itype !! Either 0, 1, 2, 3, ..., 8
        integer (ip),   intent (in)     :: isynt
        real (wp),      intent (in)     :: rsphere
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: m, n !! Counters
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy_sphere()

        ! Set constants
        this%NUMBER_OF_LATITUDES = nlat
        this%NUMBER_OF_LONGITUDES = nlon
        this%TRIANGULAR_TRUNCATION_LIMIT = ntrunc
        this%NUMBER_OF_SYNTHESES = isynt
        this%RADIUS_OF_SPHERE = rsphere

        ! Set scalar symmetries
        call this%get_scalar_symmetries( isym )

        ! Set vector symmetries
        call this%get_vector_symmetries( itype )

        ! Allocate memory
        associate( nm_dim => (ntrunc+1)*(ntrunc+2)/2 )
            allocate( this%INDEX_ORDER_M(nm_dim) )
            allocate( this%INDEX_DEGREE_N(nm_dim) )
            allocate( this%laplacian_coefficients(nm_dim) )
            allocate( this%inverse_laplacian_coefficients(nm_dim) )
            allocate( this%complex_spectral_coefficients(nm_dim) )
            allocate( this%vorticity_and_divergence_coefficients(nlat) )

            ! Fill arrays
            associate( &
                ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
                indxm => this%INDEX_ORDER_M, &
                indxn => this%INDEX_DEGREE_N, &
                lap => this%laplacian_coefficients, &
                ilap => this%inverse_laplacian_coefficients, &
                sqnn => this%vorticity_and_divergence_coefficients &
                )
                ! Set indices of order m
                indxm = [ ((m, n=m, ntrunc), m=0, ntrunc) ]
                ! Set indices of degree n
                indxn = [ ((n, n=m, ntrunc), m=0, ntrunc) ]
                ! Set laplacian coefficients
                lap = -real(indxn, kind=wp) * real(indxn + 1, kind=wp)/(rsphere**2)
                ! Set inverse laplacian coefficients
                ilap(1) = 0.0_wp
                ilap(2:nm_dim) = 1.0_wp/lap(2:nm_dim)
                ! Set vorticity and divergence coefficients
                sqnn = [ (sqrt(real((n - 1) * n, kind=wp)/rsphere), n=1, nlat) ]
            end associate
        end associate

        ! Initialize derived data types
        associate( &
            grid => this%grid, &
            workspace => this%workspace, &
            trig_func => this%trigonometric_functions, &
            unit_vectors => this%unit_vectors &
            )
            call trig_func%create( grid )
            call unit_vectors%create( grid, trig_func )
        end associate

        ! Set initialization flag
        this%initialized = .true.
        
    end subroutine create_sphere



    subroutine destroy_sphere( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if ( this%initialized .eqv. .false. ) return

         ! Release memory
        if (allocated(this%INDEX_ORDER_M)) deallocate(this%INDEX_ORDER_M)
        if (allocated(this%INDEX_DEGREE_N)) deallocate(this%INDEX_DEGREE_N)
        if (allocated(this%laplacian_coefficients )) &
            deallocate(this%laplacian_coefficients)
        if (allocated(this%inverse_laplacian_coefficients)) &
            deallocate(this%inverse_laplacian_coefficients)
        if (allocated(this%complex_spectral_coefficients)) &
            deallocate(this%complex_spectral_coefficients)
        if (allocated(this%vorticity_and_divergence_coefficients)) &
            deallocate(this%vorticity_and_divergence_coefficients)

        ! Release memory from polymorphic class variables
        if (allocated(this%grid)) deallocate(this%grid)
        if (allocated(this%workspace)) deallocate( this%workspace )

        ! Release memory from derived data types
        call this%trigonometric_functions%destroy()
        call this%unit_vectors%destroy()

        ! Reset constants
        this%NUMBER_OF_LONGITUDES = 0
        this%NUMBER_OF_LATITUDES = 0
        this%TRIANGULAR_TRUNCATION_LIMIT = 0
        this%SCALAR_SYMMETRIES = 0
        this%VECTOR_SYMMETRIES = 0
        this%NUMBER_OF_SYNTHESES = 0

        ! Reset initialization flag
        this%initialized = .false.

    end subroutine destroy_sphere



    subroutine perform_complex_analysis( this, scalar_function )
        !
        ! Purpose:
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
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out)  :: this
        real (wp),      intent (in)      :: scalar_function(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in PERFORM_COMPLEX_ANALYSIS'
        end if

        ! compute the (real) spherical harmonic coefficients
        call this%perform_scalar_analysis( scalar_function )

        ! Set complex spherical harmonic coefficients
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, & ! set the triangular truncation limit
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            psi => this%complex_spectral_coefficients &
            )
            psi = 0.5_wp * cmplx( &
                [((a(m+1, n+1), n=m, ntrunc), m=0, ntrunc)], &
                [((b(m+1, n+1), n=m, ntrunc), m=0, ntrunc)], &
                kind=wp )
        end associate

    end subroutine perform_complex_analysis



    subroutine perform_complex_synthesis( this, scalar_function )
        !
        ! Purpose:
        ! Converts gridded input array to (complex) spherical harmonic coefficients
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out)  :: this
        real (wp),      intent (out)     :: scalar_function(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m, i !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in PERFORM_COMPLEX_SYNTHESIS'
        end if

        ! Convert complex spherical harmonic coefficients to real version
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, & ! set the triangular truncation limit
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            psi => this%complex_spectral_coefficients &
            )
            ! Initialize (real) coefficients
            a = 0.0_wp
            b = 0.0_wp

            ! Fill real arrays with contents of spec
            do m = 0, ntrunc
                do n = m, ntrunc
                    ! set the spectral index
                    associate( nm => sum ([(i, i=ntrunc+1, ntrunc-m+2, -1)]) + n-m+1 )
                        ! set the real component
                        a( m+1, n+1 ) = 2.0_wp * real(psi(nm))
                        ! set the imaginary component
                        b( m+1, n+1 ) = 2.0_wp * aimag(psi(nm))
                    end associate
                end do
            end do
        end associate

        ! synthesise the scalar function from the (real) harmonic coefficients
        call this%perform_scalar_synthesis( scalar_function )

    end subroutine perform_complex_synthesis


    subroutine analyze_into_complex_spectral_coefficients( this, &
        scalar_function, spectral_coefficients )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out)  :: this
        real (wp),      intent (in)      :: scalar_function(:, :)
        complex (wp),   intent (out)     :: spectral_coefficients(:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in ANALYZE_INTO_COMPLEX_SPECTRAL_COEFFICIENTS'
        end if

        ! analyze the scalar function into (complex) harmonic coefficients
        call this%perform_complex_analysis( scalar_function )

        ! copy coefficients
        spectral_coefficients = this%complex_spectral_coefficients

    end subroutine analyze_into_complex_spectral_coefficients


    subroutine synthesize_from_complex_spectral_coefficients( this, &
        spectral_coefficients, scalar_function )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out)  :: this
        complex (wp),   intent (in)      :: spectral_coefficients(:)
        real (wp),      intent (out)     :: scalar_function(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: n, m, i !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in SYNTHESIZE_FROM_COMPLEX_SPECTRAL_COEFFICIENTS'
        end if

        ! Convert complex coefficients to real version
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, & ! set the triangular truncation limit
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            psi => spectral_coefficients &
            )
            ! Initialize (real) coefficients
            a = 0.0_wp
            b = 0.0_wp

            ! fill real arrays with contents of spec
            do m = 0, ntrunc
                do n = m, ntrunc
                     ! set the spectral index
                    associate( nm => sum ([(i, i=ntrunc+1, ntrunc-m+2, -1)]) + n-m+1 )
                        ! set the real component
                        a(m + 1, n + 1) = 2.0_wp * real(psi(nm))
                        ! set the imaginary component
                        b(m + 1, n + 1) = 2.0_wp * aimag(psi(nm))
                    end associate
                end do
            end do

        end associate

        ! synthesise the scalar function from the (real) harmonic coefficients
        call this%perform_scalar_synthesis( scalar_function )

    end subroutine synthesize_from_complex_spectral_coefficients



    subroutine perform_vector_analysis( this, vector_field )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: vector_field(:, :, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        real (wp), allocatable :: polar_component(:, :)
        real (wp), allocatable :: azimuthal_component(:, :)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): '&
                //'uninitialized object in PERFORM_VECTOR_ANALYSIS'
        end if

        ! Allocate memory
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )
            allocate( polar_component(nlat, nlon) )
            allocate( azimuthal_component(nlat, nlon) )
        end associate

        ! compute the spherical angle components
        associate( &
            F => vector_field, &
            v => polar_component, &
            w => azimuthal_component &
            )
            ! Get spherical components
            call this%unit_vectors%get_spherical_angle_components( F, v, w )
            ! Perform vector analysis
            call this%vector_analysis_from_spherical_components( v, w )
        end associate

        ! Release memory
        deallocate( polar_component)
        deallocate( azimuthal_component)

    end subroutine perform_vector_analysis



    subroutine get_scalar_laplacian( this, scalar_function, scalar_laplacian )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: scalar_function(:, :)
        real (wp),      intent (out)    :: scalar_laplacian(:, :)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_SCALAR_LAPLACIAN'
        end if

        ! Set (real) scalar spherica harmonic coefficients
        call this%perform_complex_analysis( scalar_function )

        ! Associate various quantities
        associate( &
            psi => this%complex_spectral_coefficients, &
            lap => this%laplacian_coefficients &
            )
            psi = lap * psi
        end associate

        ! Synthesize complex coefficients into gridded array
        call this%perform_complex_synthesis( scalar_laplacian )

    end subroutine get_scalar_laplacian



    subroutine invert_scalar_laplacian( this, source, solution)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: source(:, :)
        real (wp),      intent (out)    :: solution(:, :)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in INVERT_SCALAR_LAPLACIAN'
        end if

        ! Set (real) scalar spherica harmonic coefficients
        call this%perform_complex_analysis( source )

        ! Associate various quantities
        associate( &
            psi => this%complex_spectral_coefficients, &
            ilap => this%inverse_laplacian_coefficients &
            )
            psi = ilap * psi
        end associate

        ! Synthesize complex coefficients into gridded array
        call this%perform_complex_synthesis( solution )

    end subroutine invert_scalar_laplacian


    subroutine invert_helmholtz( this, helmholtz_constant, source, solution )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), target, intent (in out) :: this
        real (wp),              intent (in)     :: helmholtz_constant
        real (wp),              intent (in)     :: source(:, :)
        real (wp),              intent (out)    :: solution(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        real (wp), parameter :: ZERO=nearest(1.0_wp, 1.0_wp)-nearest(1.0_wp, -1.0_wp)
        real (wp), pointer   :: iptr(:) => null()
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in INVERT_HELMHOLTZ'
        end if

        ! Set (real) scalar spherica harmonic coefficients
        call this%perform_complex_analysis( source )

        ! Associate various quantities
        associate( &
            nm_dim => size(this%complex_spectral_coefficients), &
            psi => this%complex_spectral_coefficients, &
            ilap => this%inverse_laplacian_coefficients, &
            lap => this%laplacian_coefficients, &
            xlmbda => helmholtz_constant &
            )
            if ( xlmbda == ZERO ) then
                ! Associate pointer
                iptr => ilap
            else
                ! Allocate memory
                allocate( iptr(nm_dim) )
                iptr = 1.0_wp/(lap - xlmbda)
            end if

            ! Set coefficients
            psi = iptr * psi

            ! Synthesize complex coefficients into gridded array
            call this%perform_complex_synthesis( solution )

            ! Garbage collection
            if ( xlmbda == ZERO ) then
                ! Nullify pointer
                nullify( iptr )
            else
                ! Release memory
                deallocate( iptr )
            end if
        end associate

    end subroutine invert_helmholtz



    subroutine get_gradient( this, scalar_function, &
        polar_gradient_component, azimuthal_gradient_component )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: scalar_function(:, :)
        real (wp),      intent (out)    :: polar_gradient_component(:, :)
        real (wp),      intent (out)    :: azimuthal_gradient_component(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_GRADIENT'
        end if

        ! compute the (real) spherical harmonic coefficients
        call this%perform_scalar_analysis( scalar_function )

        ! compute gradient
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            v => polar_gradient_component, &
            w => azimuthal_gradient_component, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            ! Initialize polar coefficients
            br = 0.0_wp
            bi = 0.0_wp

            ! Set polar coefficients
            do n=1, nlat
                br(:, n) = a(:, n) * sqnn(n)
                bi(:, n) = b(:, n) * sqnn(n)
            end do
            ! Set azimuthal coefficients
            cr = 0.0_wp
            ci = 0.0_wp
            ! Compute vector harmonic synthesis
            call this%perform_vector_synthesis( v, w )
        end associate

    end subroutine get_gradient


    subroutine invert_gradient( this, source, solution )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: source(:, :, :)
        real (wp),      intent (out)    :: solution(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in INVERT_GRADIENT'
        end if

        ! compute the (real) spherical harmonic coefficients
        call this%perform_vector_analysis( source )

        ! Invert gradient
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            f => solution, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients &
            )
            ! Initialize (real) coefficients
            a = 0.0_wp
            b = 0.0_wp

            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    ! Set real coefficients
                    a(m, n) = sqnn(n) * br(m, n)
                    ! Set imaginary coefficients
                    b(m, n) = sqnn(n) * bi(m, n)
                end do
            end do
            ! Compute scalar synthesis
            call this%perform_scalar_synthesis( f )
        end associate

    end subroutine invert_gradient


    subroutine get_vorticity( this, vector_field, vorticity )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: vector_field(:, :, :)
        real (wp),      intent (out)    :: vorticity(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_VORTICITY'
        end if

        ! Calculate the (real) vector harmonic coefficients
        call this%perform_vector_analysis( vector_field )

        ! compute gradient
        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            vort => vorticity, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            ! Initialize (real) coefficients
            a = 0.0_wp
            b = 0.0_wp

            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    ! Set real coefficients
                    a(m, n) = sqnn(n) * cr(m, n)
                    ! Set imaginary coefficients
                    b(m, n) = sqnn(n) * ci(m, n)
                end do
            end do
            ! Compute vector harmonic synthesis
            call this%perform_scalar_synthesis( vort )
        end associate


    end subroutine get_vorticity


    subroutine invert_vorticity( this, source, solution )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: source(:, :, :)
        real (wp),      intent (out)    :: solution(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in INVERT_VORTICITY'
        end if

        ! compute the (real) spherical harmonic coefficients
        call this%perform_vector_analysis( source )

        ! Invert gradient
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            f => solution, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            ! Initialize azimuthal coefficients
            cr = 0.0_wp
            ci = 0.0_wp

            ! Compute m = 0 coefficients
            do n = 1, nlat
                ! Set real coefficients
                cr(1, n) = a(1, n)/sqnn(n)
                ! Set imaginary coeffients
                ci(1, n) = b(1, n)/sqnn(n)
            end do

            ! Compute m >0 coefficients
            do m=2, ntrunc+1
                do n=m, ntrunc+1
                    ! Set real coefficients
                    cr(m, n) = a(m, n)/sqnn(n)
                    ! Set imaginary coefficients
                    ci(m, n) = b(m, n)/sqnn(n)
                end do
            end do
            ! Compute scalar synthesis
            call this%perform_scalar_synthesis( f )
        end associate

    end subroutine invert_vorticity

    subroutine get_divergence_from_coefficients(this, divergence)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (out)    :: divergence (:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_DIVERGENCE_FROM_COEFFICIENTS'
        end if

        ! compute gradient
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            div => divergence, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients &
            )
            ! Initialize (real) coefficients
            a = 0.0_wp
            b = 0.0_wp

            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    ! Set real coefficients
                    a(m, n) = -sqnn(n) * br(m, n)
                    ! Set imaginary coefficients
                    b(m, n) = -sqnn(n) * bi(m, n)
                end do
            end do
            ! Compute vector harmonic synthesis
            call this%perform_scalar_synthesis( div )
        end associate

    end subroutine get_divergence_from_coefficients



    subroutine get_divergence_from_vector_field( this, vector_field, divergence )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: vector_field (:, :, :)
        real (wp),      intent (out)    :: divergence (:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_DIVERGENCE_FROM_VECTOR_FIELD'
        end if

        ! Calculate the (real) vector harmonic coefficients
        call this%perform_vector_analysis( vector_field )

        ! compute divergence
        call get_divergence_from_coefficients(this, divergence)

    end subroutine get_divergence_from_vector_field


    subroutine get_divergence_from_spherical_angles( this, &
        polar_component, azimuthal_component, divergence )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: polar_component(:,:)
        real (wp),      intent (in)     :: azimuthal_component(:,:)
        real (wp),      intent (out)    :: divergence (:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_DIVERGENCE_FROM_SPHERICAL_ANGLES'
        end if

        ! Calculate the (real) vector harmonic coefficients
        associate( &
            v => polar_component, &
            w => azimuthal_component &
            )
            call this%vector_analysis_from_spherical_components( v, w )
        end associate

        ! compute divergence
        call get_divergence_from_coefficients(this, divergence)

    end subroutine get_divergence_from_spherical_angles



    subroutine invert_divergence(this, source, polar_solution, azimuthal_solution )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: source(:,:)
        real (wp),      intent (out)    :: polar_solution(:,:)
        real (wp),      intent (out)    :: azimuthal_solution(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip):: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in INVERT_DIVERGENCE'
        end if

        ! Calculate the (real) scalar harmonic coefficients
        associate( div => source )
            call this%perform_scalar_analysis( div )
        end associate

        ! Invert gradient
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            v => polar_solution, &
            w => azimuthal_solution, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            ! Initialize polar coefficients
            br = 0.0_wp
            bi = 0.0_wp

            ! Compute m = 0 coefficients
            do n = 1, nlat
                ! Set real coefficients
                br(1, n) = -a(1, n)/sqnn(n)
                ! Set imaginary coeffients
                bi(1, n) = -b(1, n)/sqnn(n)
            end do

            ! Compute m >0 coefficients
            do m=2, ntrunc+1
                do n=m, ntrunc+1
                    ! Set real coefficients
                    br(m, n) = -a(m, n)/sqnn(n)
                    ! Set imaginary coefficients
                    bi(m, n) = -b(m, n)/sqnn(n)
                end do
            end do
            ! Set azimuthal coefficients
            cr = 0.0_wp !Not present in original source
            ci = 0.0_wp
            ! Compute scalar synthesis
            call this%perform_vector_synthesis(v, w)
        end associate

    end subroutine invert_divergence


    subroutine get_vorticity_and_divergence_coefficients_from_velocities(this, &
        polar_component, azimuthal_component, vort_spec, div_spec)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere),         intent (in out) :: this
        real (wp),              intent (in)     :: polar_component(:, :)
        real (wp),              intent (in)     :: azimuthal_component(:, :)
        complex (wp),           intent (out)    :: vort_spec(:)
        complex (wp),           intent (out)    :: div_spec(:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: n, m  !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_VORTICITY_AND_DIVERGENCE_COEFFICIENTS_FROM_VELOCITIES'
        end if

        ! Perform vector analysis
        call this%vector_analysis_from_spherical_components( &
            polar_component, azimuthal_component )

        ! Associate various quantities
        associate( &
            nm_dim => size( this%complex_spectral_coefficients ), &
            nlat => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            rsphere => this%RADIUS_OF_SPHERE, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            ! Initialize real scalar coefficients
            a = 0.0_wp
            b = 0.0_wp
            ! Set auxiliary coefficients for divergence
            do n = 1, nlat
                a(:, n) = -sqnn(n) * br(:, n)
                b(:, n) = -sqnn(n) * bi(:, n)
            end do
            ! Set divergence coefficients
            div_spec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp )
            ! Set auxiliary coefficients for vorticity
            do n=1, nlat
                a(:, n) = sqnn(n) * cr(:, n)
                b(:, n) = sqnn(n) * ci(:, n)
            end do
            ! Set vorticity coefficients
            vort_spec = 0.5_wp * cmplx( &
                [((a(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                [((b(m, n), n=m, ntrunc+1), m=1, ntrunc+1)], &
                kind=wp )
        end associate

    end subroutine get_vorticity_and_divergence_coefficients_from_velocities


    subroutine get_velocities_from_vorticity_and_divergence_coefficients(this, &
        vort_spec, div_spec, polar_component, azimuthal_component)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        complex (wp),   intent (in)     :: vort_spec(:)
        complex (wp),   intent (in)     :: div_spec(:)
        real (wp),      intent (out)    :: polar_component(:, :)
        real (wp),      intent (out)    :: azimuthal_component(:, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: nm, n, m, i  !! Counters
        real (wp), allocatable :: isqnn(:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in GET_VELOCITIES_FROM_VORTICITY_AND_DIVERGENCE_COEFFICIENTS'
        end if

        ! Associate various quantities
        associate( &
            v => polar_component, &
            w => azimuthal_component, &
            nlat => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            sqnn => this%vorticity_and_divergence_coefficients, &
            a => this%workspace%real_harmonic_coefficients, &
            b => this%workspace%imaginary_harmonic_coefficients, &
            br => this%workspace%real_polar_harmonic_coefficients, &
            bi => this%workspace%imaginary_polar_harmonic_coefficients, &
            cr => this%workspace%real_azimuthal_harmonic_coefficients, &
            ci => this%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            ! Allocate memory
            allocate( isqnn(size(sqnn)) )

            isqnn(1) = 0.0_wp
            do n = 2, nlat
                isqnn(n) = 1.0_wp/sqnn(n)
            end do

            ! Initialize (real) coefficients
            a = 0.0_wp
            b = 0.0_wp

            ! Set (real) coefficients
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = -2.0_wp * real(div_spec(nm))
                    b(m, n) = -2.0_wp * aimag(div_spec(nm))
                end do
            end do

            ! Initialize polar coefficients
            br = 0.0_wp
            bi = 0.0_wp

            ! Set polar coefficients
            do n=1, nlat
                br(:, n) = isqnn(n)*a(:, n)
                bi(:, n) = isqnn(n)*b(:, n)
            end do

            ! Re-initialize (real) coefficients
            a = 0.0_wp
            b = 0.0_wp

            ! Set azimuthal coefficients
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = 2.0_wp * real(vort_spec(nm))
                    b(m, n) = 2.0_wp * aimag(vort_spec(nm))
                end do
            end do

            ! Initialize azimuthal coefficients
            cr = 0.0_wp
            ci = 0.0_wp

            ! Set azimuthal coefficients
            do n=1, nlat
                cr(:, n) = isqnn(n)*a(:, n)
                ci(:, n) = isqnn(n)*b(:, n)
            end do

            ! Compute vector harmonic sysnthesis to get components
            call this%perform_vector_synthesis( v, w )
        end associate

        ! Release memory
        deallocate( isqnn )

    end subroutine get_velocities_from_vorticity_and_divergence_coefficients



    subroutine compute_angular_momentum(this, scalar_function, angular_momentum)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        real (wp),      intent (in)     :: scalar_function(:, :)
        real (wp),      intent (out)    :: angular_momentum(:, :, :)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: k, l   !! Counters
        real (wp), allocatable :: polar_gradient_component(:, :)
        real (wp), allocatable :: azimuthal_gradient_component(:, :)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Sphere): uninitialized object'&
                //' in COMPUTE_ANGULAR_MOMENTUM'
        end if

        ! Allocate memory
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES &
            )
            allocate( polar_gradient_component(nlat, nlon) )
            allocate( azimuthal_gradient_component(nlat, nlon) )
        end associate

        ! Calculate the spherical surface gradient components
        associate( &
            f => scalar_function, &
            grad_theta => polar_gradient_component, &
            grad_phi => azimuthal_gradient_component &
            )
            call this%get_gradient( f, grad_theta, grad_phi )
        end associate

        ! Calculate the rotation operator applied to a scalar function
        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES, &
            R => angular_momentum &
            )
            do l = 1, nlon
                do k = 1, nlat
                    associate( &
                        theta => this%unit_vectors%polar(k, l), &
                        phi => this%unit_vectors%azimuthal( k, l), &
                        grad_theta => polar_gradient_component(k, l), &
                        grad_phi => azimuthal_gradient_component(k, l) &
                        )
                        R(:, k, l) = phi * grad_theta - theta * grad_phi
                    end associate
                end do
            end do
        end associate

        ! Release memory
        deallocate( polar_gradient_component )
        deallocate( azimuthal_gradient_component)

    end subroutine compute_angular_momentum
    


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
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out)  :: this
        integer (ip),   intent (in)      :: n
        integer (ip),   intent (in)      :: m
        integer (ip)                     :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: i !! Counter
        !----------------------------------------------------------------------

        ! Initialize return value
        return_value = -1

        associate( ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT )
            if ( m <= n .and. max(n, m) <= ntrunc ) then
                return_value = sum ([ (i, i=ntrunc+1, ntrunc-m+2, -1) ]) + n-m+1
            end if
        end associate

    end function get_index


    function get_coefficient( this, n, m ) result ( return_value )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        integer (ip),   intent (in)     :: n
        integer (ip),   intent (in)     :: m
        complex (wp)                    :: return_value
        !----------------------------------------------------------------------

        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            nm  => this%get_index(n, m), &
            nm_conjg => this%get_index(n, -m), &
            psi => this%complex_spectral_coefficients &
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


    subroutine get_scalar_symmetries( this, isym )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        integer (ip),   intent (in)     :: isym
        !----------------------------------------------------------------------

        select case (isym)
            case (2)
                this%SCALAR_SYMMETRIES = isym
            case (1)
                this%SCALAR_SYMMETRIES = isym
            case (0)
                this%SCALAR_SYMMETRIES = isym
            case default
                error stop 'TYPE (Sphere) in get_SCALAR_SYMMETRIES'&
                    //'Optional argument isym must be either '&
                    //' 0, 1, or 2 (default isym = 0)'
        end select

    end subroutine get_scalar_symmetries
    

    subroutine get_vector_symmetries( this, itype )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Sphere), intent (in out) :: this
        integer (ip),   intent (in)     :: itype
        !----------------------------------------------------------------------

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
                error stop 'TYPE (Sphere) in GET_VECTOR_SYMMETRIES'&
                    //'Optional argument itype must be either '&
                    //' 0, 1, 2, ..., 8 (default itype = 0)'
        end select

    end subroutine get_vector_symmetries


end module type_Sphere
