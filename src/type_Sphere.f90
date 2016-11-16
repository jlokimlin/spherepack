module type_Sphere

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Workspace, only: &
        Workspace

    use type_SphericalGrid, only: &
        SphericalGrid

    use type_SphericalUnitVectors, only: &
        SphericalUnitVectors

    use type_TrigonometricFunctions, only: &
        TrigonometricFunctions

    use type_Vector3D, only: &
        Vector => Vector3D, &
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
        ! Type components
        !----------------------------------------------------------------------
        logical,                            public :: initialized = .false.
        integer(ip),                        public :: NUMBER_OF_LONGITUDES = 0
        integer(ip),                        public :: NUMBER_OF_LATITUDES = 0
        integer(ip),                        public :: TRIANGULAR_TRUNCATION_LIMIT = 0
        integer(ip),                        public :: SCALAR_SYMMETRIES = 0
        integer(ip),                        public :: VECTOR_SYMMETRIES = 0
        integer(ip),                        public :: NUMBER_OF_SYNTHESES = 0
        integer(ip),          allocatable,  public :: INDEX_ORDER_M(:)
        integer(ip),          allocatable,  public :: INDEX_DEGREE_N(:)
        real(wp),                           public :: RADIUS_OF_SPHERE = 0.0_wp
        real(wp),             allocatable,  public :: vorticity_and_divergence_coefficients(:)
        real(wp),             allocatable,  public :: laplacian_coefficients(:)
        real(wp),             allocatable,  public :: inverse_laplacian_coefficients(:)
        complex(wp),          allocatable, public :: complex_spectral_coefficients(:)
        class(Workspace),     allocatable, public :: workspace
        class(SphericalGrid), allocatable, public :: grid
        type(TrigonometricFunctions),      public :: trigonometric_functions
        type(SphericalUnitVectors),        public :: unit_vectors
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
        !----------------------------------------------------------------------
        procedure, public  :: create_sphere
        procedure, public  :: destroy => destroy_sphere
        procedure, public  :: destroy_sphere
        procedure, public  :: get_index
        procedure, public  :: get_coefficient
        procedure, public  :: invert_helmholtz
        procedure, public  :: invert_vorticity
        procedure, public  :: get_gradient
        procedure, private :: invert_gradient_from_spherical_components
        procedure, private :: get_vorticity_from_spherical_components
        procedure, private :: get_vorticity_from_vector_field
        procedure, private :: get_divergence_from_vector_field
        procedure, private :: get_divergence_from_spherical_components
        procedure, public  :: invert_divergence
        procedure, public  :: get_rotation_operator => compute_angular_momentum
        procedure, private :: get_scalar_symmetries
        procedure, private :: get_vector_symmetries
        procedure, public  :: perform_complex_analysis
        procedure, public  :: perform_complex_synthesis
        procedure, private :: perform_vector_analysis_from_vector_field
        procedure, public  :: synthesize_from_complex_spectral_coefficients
        procedure, public  :: analyze_into_complex_spectral_coefficients
        procedure, public  :: get_vorticity_and_divergence_from_velocities
        procedure, private :: get_velocities_from_vorticity_and_divergence_coefficients
        procedure, public  :: get_velocities_from_vorticity_and_divergence
        procedure, private :: get_scalar_laplacian
        procedure, private :: compute_vector_laplacian_coefficients
        procedure, private :: get_vector_laplacian_from_spherical_components
        procedure, private :: get_vector_laplacian_from_vector_field
        procedure, private :: invert_scalar_laplacian
        procedure, private :: invert_vector_laplacian
        !----------------------------------------------------------------------
        ! Deferred type-bound procedures
        !----------------------------------------------------------------------
        procedure (scalar_analysis),  deferred, public :: &
            perform_scalar_analysis
        procedure (scalar_synthesis), deferred, public :: &
            perform_scalar_synthesis
        procedure (vector_analysis),  deferred, public :: &
            vector_analysis_from_spherical_components
        procedure (vector_synthesis), deferred, public :: &
            perform_vector_synthesis
        !----------------------------------------------------------------------
        ! Generic type-bound procedures
        !----------------------------------------------------------------------
        generic, public :: perform_vector_analysis => &
            perform_vector_analysis_from_vector_field
        generic, public :: invert_laplacian => &
            invert_scalar_laplacian, &
            invert_vector_laplacian
        generic, public :: get_divergence => &
            get_divergence_from_vector_field, &
            get_divergence_from_spherical_components
        generic, public :: get_laplacian => &
            get_scalar_laplacian, &
            get_vector_laplacian_from_spherical_components, &
            get_vector_laplacian_from_vector_field
        generic, public :: get_vorticity => &
            get_vorticity_from_spherical_components, &
            get_vorticity_from_vector_field
        generic, public :: invert_gradient => &
            invert_gradient_from_spherical_components
        !----------------------------------------------------------------------
    end type Sphere

    abstract interface
        subroutine scalar_analysis(self, scalar_function)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(in)     :: scalar_function(:,:)
            !----------------------------------------------------------------------
        end subroutine scalar_analysis

        subroutine scalar_synthesis(self, scalar_function)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(out)   :: scalar_function(:,:)
            !----------------------------------------------------------------------
        end subroutine scalar_synthesis

        subroutine vector_analysis(self, polar_component, azimuthal_component)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(in)     :: polar_component(:,:)
            real(wp),      intent(in)     :: azimuthal_component(:,:)
            !----------------------------------------------------------------------
        end subroutine vector_analysis

        subroutine vector_synthesis(self, polar_component, azimuthal_component)
            import :: Sphere, wp
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            class(Sphere), intent(inout)  :: self
            real(wp),      intent(out)    :: polar_component(:,:)
            real(wp),      intent(out)    :: azimuthal_component(:,:)
            !----------------------------------------------------------------------
        end subroutine vector_synthesis
    end interface

contains

    subroutine create_sphere(self, nlat, nlon, ntrunc, isym, itype, nt, rsphere)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        integer(ip),   intent(in)     :: nlat
        integer(ip),   intent(in)     :: nlon
        integer(ip),   intent(in)     :: ntrunc
        integer(ip),   intent(in)     :: isym  !! Either 0, 1, or 2
        integer(ip),   intent(in)     :: itype !! Either 0, 1, 2, 3, ..., 8
        integer(ip),   intent(in)     :: nt
        real(wp),      intent(in)     :: rsphere
        !--------------------------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------------------------
        integer(ip) :: m, n !! Counters
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call self%destroy_sphere()

        !
        !==> Set constants
        !
        self%NUMBER_OF_LATITUDES = nlat
        self%NUMBER_OF_LONGITUDES = nlon
        self%TRIANGULAR_TRUNCATION_LIMIT = ntrunc
        self%NUMBER_OF_SYNTHESES = nt
        self%RADIUS_OF_SPHERE = rsphere

        !
        !==> Set scalar symmetries
        !
        call self%get_scalar_symmetries(isym)

        !
        !==> Set vector symmetries
        !
        call self%get_vector_symmetries(itype)

        !
        !==> Allocate memory
        !
        associate( nm_dim => (ntrunc+1)*(ntrunc+2)/2 )
            allocate(self%INDEX_ORDER_M(nm_dim) )
            allocate(self%INDEX_DEGREE_N(nm_dim) )
            allocate(self%laplacian_coefficients(nm_dim) )
            allocate(self%inverse_laplacian_coefficients(nm_dim) )
            allocate(self%complex_spectral_coefficients(nm_dim) )
            allocate(self%vorticity_and_divergence_coefficients(nlat) )
            !
            !==> Fill arrays
            !
            associate( &
                ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
                indxm => self%INDEX_ORDER_M, &
                indxn => self%INDEX_DEGREE_N, &
                lap => self%laplacian_coefficients, &
                ilap => self%inverse_laplacian_coefficients, &
                sqnn => self%vorticity_and_divergence_coefficients &
                )
                !
                !==> Precompute indices of order m
                !
                indxm = [ ((m, n=m, ntrunc), m=0, ntrunc) ]
                !
                !==> Precompute indices of degree n
                !
                indxn = [ ((n, n=m, ntrunc), m=0, ntrunc) ]
                !
                !==> Precompute laplacian coefficients
                !
                lap = -real(indxn, kind=wp) * real(indxn + 1, kind=wp)/(rsphere**2)
                !
                !==> Precompute inverse laplacian coefficients
                !
                ilap(1) = 0.0_wp
                ilap(2:nm_dim) = 1.0_wp/lap(2:nm_dim)
                !
                !==> Precompute vorticity and divergence coefficients
                !
                sqnn = [ (sqrt(real((n - 1) * n, kind=wp)/rsphere), n=1, nlat) ]

            end associate
        end associate
        !
        !==> Initialize derived data types
        !
        associate( &
            grid => self%grid, &
            trig_func => self%trigonometric_functions, &
            unit_vectors => self%unit_vectors &
            )

            trig_func = TrigonometricFunctions(grid)

            unit_vectors = SphericalUnitVectors(grid)

        end associate

        !
        !==> Set initialization flag
        !
        self%initialized = .true.
        
    end subroutine create_sphere



    subroutine destroy_sphere(self)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        !----------------------------------------------------------------------

        ! Check flag
        if (.not.self%initialized) return

        !
        !==> Release memory
        !
        if (allocated(self%INDEX_ORDER_M)) then
            deallocate(self%INDEX_ORDER_M)
        end if

        if (allocated(self%INDEX_DEGREE_N)) then
            deallocate(self%INDEX_DEGREE_N)
        end if

        if (allocated(self%laplacian_coefficients )) then
            deallocate(self%laplacian_coefficients)
        end if

        if (allocated(self%inverse_laplacian_coefficients)) then
            deallocate(self%inverse_laplacian_coefficients)
        end if

        if (allocated(self%complex_spectral_coefficients)) then
            deallocate(self%complex_spectral_coefficients)
        end if

        if (allocated(self%vorticity_and_divergence_coefficients)) then
            deallocate(self%vorticity_and_divergence_coefficients)
        end if

        !
        !==> Release memory from polymorphic class variables
        !
        if (allocated(self%grid)) deallocate( self%grid )

        if (allocated(self%workspace)) deallocate( self%workspace )
        !
        !==>  Release memory from derived data types
        !
        call self%trigonometric_functions%destroy()
        call self%unit_vectors%destroy()

        !
        !==> Reset constants
        !
        self%NUMBER_OF_LONGITUDES = 0
        self%NUMBER_OF_LATITUDES = 0
        self%TRIANGULAR_TRUNCATION_LIMIT = 0
        self%SCALAR_SYMMETRIES = 0
        self%VECTOR_SYMMETRIES = 0
        self%NUMBER_OF_SYNTHESES = 0

        !
        !==> Reset initialization flag
        !
        self%initialized = .false.

    end subroutine destroy_sphere



    subroutine perform_complex_analysis(self, scalar_function)
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
        ! In Fortran syntax, values of m (degree) and n (order) as a function
        ! of the index nm are:
        !
        ! integer(ip), dimension ((mtrunc+1)*(mtrunc+2)/2) :: indxm, indxn
        ! indxm = [((m, n=m, mtrunc), m=0, mtrunc)]
        ! indxn = [((n, n=m, mtrunc), m=0, mtrunc)]
        !
        ! Conversely, the index nm as a function of m and n is:
        ! nm = sum([(i, i=mtrunc+1, mtrunc-m+2, -1)])+n-m+1
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in perform_complex_analysis'
        end if

        !
        !==>  compute the (real) spherical harmonic coefficients
        !
        associate( f => scalar_function )

            call self%perform_scalar_analysis(f)

        end associate

        !
        !==> Set complex spherical harmonic coefficients
        !
        associate( &
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            psi => self%complex_spectral_coefficients &
            )

            psi = 0.5_wp * cmplx( &
                [((a(m+1, n+1), n=m, ntrunc), m=0, ntrunc)], &
                [((b(m+1, n+1), n=m, ntrunc), m=0, ntrunc)], &
                kind=wp )

        end associate

    end subroutine perform_complex_analysis



    subroutine perform_complex_synthesis(self, scalar_function )
        !
        ! Purpose:
        ! Converts gridded input array to (complex) spherical harmonic coefficients
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(out)    :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m, i !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in perform_complex_synthesis'
        end if

        !
        !==> Convert complex spherical harmonic coefficients to real version
        !
        associate( &
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            psi => self%complex_spectral_coefficients &
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

        !
        !==> synthesise the scalar function from the (real) harmonic coefficients
        !
        call self%perform_scalar_synthesis(scalar_function)

    end subroutine perform_complex_synthesis



    subroutine analyze_into_complex_spectral_coefficients(self, &
        scalar_function, spectral_coefficients )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: scalar_function(:,:)
        complex(wp),   intent(out)    :: spectral_coefficients(:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in analyze_into_complex_spectral_coefficients'
        end if

        ! analyze the scalar function into (complex) harmonic coefficients
        call self%perform_complex_analysis(scalar_function)

        ! copy coefficients
        spectral_coefficients = self%complex_spectral_coefficients

    end subroutine analyze_into_complex_spectral_coefficients



    subroutine synthesize_from_complex_spectral_coefficients(self, &
        spectral_coefficients, scalar_function )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        complex(wp),   intent(in)     :: spectral_coefficients(:)
        real(wp),      intent(out)    :: scalar_function(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m, i !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in synthesize_from_complex_spectral_coefficients'
        end if

        ! Convert complex coefficients to real version
        associate( &
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, & ! set the triangular truncation limit
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
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
        call self%perform_scalar_synthesis(scalar_function)

    end subroutine synthesize_from_complex_spectral_coefficients



    subroutine perform_vector_analysis_from_vector_field(self, vector_field )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)    :: vector_field(:,:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)           :: error_flag
        real(wp), allocatable :: polar_component(:,:)
        real(wp), allocatable :: azimuthal_component(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //'in perform_vector_analysis'
        end if

        ! Allocate memory
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES &
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
            call self%unit_vectors%get_spherical_angle_components(F, v, w)
            ! Perform vector analysis
            call self%vector_analysis_from_spherical_components(v, w)
        end associate

        ! Release memory
        deallocate( polar_component)
        deallocate( azimuthal_component)

    end subroutine perform_vector_analysis_from_vector_field



    subroutine get_scalar_laplacian(self, scalar_function, scalar_laplacian )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: scalar_function(:,:)
        real(wp),      intent(out)    :: scalar_laplacian(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_scalar_laplacian'
        end if

        ! Set (real) scalar spherica harmonic coefficients
        call self%perform_complex_analysis(scalar_function)

        ! Associate various quantities
        associate( &
            psi => self%complex_spectral_coefficients, &
            lap => self%laplacian_coefficients &
            )
            psi = lap * psi
        end associate

        ! Synthesize complex coefficients into gridded array
        call self%perform_complex_synthesis(scalar_laplacian)

    end subroutine get_scalar_laplacian



    subroutine invert_scalar_laplacian(self, source, solution)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)    :: source(:,:)
        real(wp),      intent(out)   :: solution(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in invert_scalar_laplacian'
        end if

        ! Set (real) scalar spherica harmonic coefficients
        call self%perform_complex_analysis( source )

        ! Associate various quantities
        associate( &
            psi => self%complex_spectral_coefficients, &
            ilap => self%inverse_laplacian_coefficients &
            )
            psi = ilap * psi
        end associate

        ! Synthesize complex coefficients into gridded array
        call self%perform_complex_synthesis( solution )

    end subroutine invert_scalar_laplacian

    subroutine compute_vector_laplacian_coefficients(self)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip) :: n !! Counter
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in compute_vector_laplacian_coefficients'
        end if

        ! compute vector laplacian
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            ityp => self%VECTOR_SYMMETRIES, &
            lap => self%laplacian_coefficients, &
            br => self%workspace%real_polar_harmonic_coefficients, &
            bi => self%workspace%imaginary_polar_harmonic_coefficients, &
            cr => self%workspace%real_azimuthal_harmonic_coefficients, &
            ci => self%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            select case (ityp)
                case (0, 3, 6)
                    !
                    !==>  All coefficients needed
                    !
                    do n=1, nlat
                        ! Set polar coefficients
                        br(:,n) = lap(n) * br(:,n)
                        bi(:,n) = lap(n) * bi(:,n)
                        ! Set azimuthal coefficients
                        cr(:,n) = lap(n) * cr(:,n)
                        ci(:,n) = lap(n) * ci(:,n)
                    end do
                case (1, 4, 7)
                    !
                    !==>     vorticity is zero so cr,ci=0 not used
                    !
                    do n=1, nlat
                        ! Set polar coefficients
                        br(:,n) = lap(n) * br(:,n)
                        bi(:,n) = lap(n) * bi(:,n)
                    end do
                case default
                    !
                    !==> divergence is zero so br,bi=0 not used
                    !
                    do n=1, nlat
                        ! Set azimuthal coefficients
                        cr(:,n) = lap(n) * cr(:,n)
                        ci(:,n) = lap(n) * ci(:,n)
                    end do
            end select
        end associate

    end subroutine compute_vector_laplacian_coefficients

    subroutine get_vector_laplacian_from_spherical_components(self, &
        polar_component, azimuthal_component, polar_laplacian, azimuthal_laplacian )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: polar_component(:,:)
        real(wp),      intent(in)     :: azimuthal_component(:,:)
        real(wp),      intent(out)    :: polar_laplacian(:,:)
        real(wp),      intent(out)    :: azimuthal_laplacian(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_vector_laplacian_from_spherical_components'
        end if

        ! Set vector spherical harmonic coefficients
        associate( &
            v => polar_component, &
            w => azimuthal_component &
            )
            call self%vector_analysis_from_spherical_components(v, w)
        end associate

        ! Compute vector laplacian coefficients
        call self%compute_vector_laplacian_coefficients()


        ! Synthesize vector laplacian from coefficients
        associate( &
            vlap => polar_laplacian, &
            wlap => azimuthal_laplacian &
            )
            call self%perform_vector_synthesis(vlap, wlap)
        end associate

    end subroutine get_vector_laplacian_from_spherical_components


    subroutine get_vector_laplacian_from_vector_field(self, &
        vector_field, polar_laplacian, azimuthal_laplacian)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: vector_field(:,:,:)
        real(wp),      intent(out)    :: polar_laplacian(:,:)
        real(wp),      intent(out)    :: azimuthal_laplacian(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_vector_laplacian_from_vector_field'
        end if

        ! Compute vector laplacian coefficients
        call self%perform_vector_analysis( vector_field )


        ! Synthesize vector laplacian from coefficients
        associate( &
            vlap => polar_laplacian, &
            wlap => azimuthal_laplacian &
            )
            call self%perform_vector_synthesis( vlap, wlap )
        end associate

    end subroutine get_vector_laplacian_from_vector_field


    subroutine invert_vector_laplacian(self, &
        polar_source, azimuthal_source, polar_solution, azimuthal_solution)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: polar_source(:,:)
        real(wp),      intent(in)     :: azimuthal_source(:,:)
        real(wp),      intent(out)    :: polar_solution(:,:)
        real(wp),      intent(out)    :: azimuthal_solution(:,:)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip) :: n !! Counter
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in invert_vector_laplacian'
        end if

        ! Set vector spherical harmonic coefficients
        associate( &
            vlap => polar_source, &
            wlap => azimuthal_source &
            )
            call self%vector_analysis_from_spherical_components( vlap, wlap )
        end associate

        ! compute vector laplacian
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            ityp => self%VECTOR_SYMMETRIES, &
            ilap => self%inverse_laplacian_coefficients, &
            br => self%workspace%real_polar_harmonic_coefficients, &
            bi => self%workspace%imaginary_polar_harmonic_coefficients, &
            cr => self%workspace%real_azimuthal_harmonic_coefficients, &
            ci => self%workspace%imaginary_azimuthal_harmonic_coefficients, &
            v => polar_solution, &
            w => azimuthal_solution &
            )
            select case (ityp)
                case (0, 3, 6)
                    !
                    !==>  All coefficients needed
                    !
                    do n=1, nlat
                        ! Set polar coefficients
                        br(:,n) = ilap(n) * br(:,n)
                        bi(:,n) = ilap(n) * bi(:,n)
                        ! Set azimuthal coefficients
                        cr(:,n) = ilap(n) * cr(:,n)
                        ci(:,n) = ilap(n) * ci(:,n)
                    end do
                case (1, 4, 7)
                    !
                    !==>     vorticity is zero so cr,ci=0 not used
                    !
                    do n=1, nlat
                        ! Set polar coefficients
                        br(:,n) = ilap(n) * br(:,n)
                        bi(:,n) = ilap(n) * bi(:,n)
                    end do
                case default
                    !
                    !==> divergence is zero so br,bi=0 not used
                    !
                    do n=1, nlat
                        ! Set azimuthal coefficients
                        cr(:,n) = ilap(n) * cr(:,n)
                        ci(:,n) = ilap(n) * ci(:,n)
                    end do
            end select
            !
            !==> synthesize coefficients inot vector field (v,w)
            !
            call self%perform_vector_synthesis(v, w)
        end associate


    end subroutine invert_vector_laplacian


    subroutine invert_helmholtz(self, helmholtz_constant, source, solution )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), target, intent(inout)  :: self
        real(wp),              intent(in)     :: helmholtz_constant
        real(wp),              intent(in)     :: source(:,:)
        real(wp),              intent(out)    :: solution(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        real(wp), parameter :: ZERO=nearest(1.0_wp, 1.0_wp)-nearest(1.0_wp, -1.0_wp)
        real(wp), pointer   :: iptr(:) => null()
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in invert_helmholtz'
        end if

        !
        !==> Set (real) scalar spherica harmonic coefficients
        !
        call self%perform_complex_analysis(source)

        ! Associate various quantities
        associate( &
            nm_dim => size(self%complex_spectral_coefficients), &
            psi => self%complex_spectral_coefficients, &
            ilap => self%inverse_laplacian_coefficients, &
            lap => self%laplacian_coefficients, &
            xlmbda => helmholtz_constant &
            )
            !
            !==> Associate local pointer
            !
            if ( xlmbda == ZERO ) then
                ! Assign pointer
                iptr => ilap
            else
                ! Allocate memory
                allocate( iptr(nm_dim) )
                iptr = 1.0_wp/(lap - xlmbda)
            end if

            !
            !==> Set coefficients
            !
            psi = iptr * psi

            !
            !==> Synthesize complex coefficients into gridded array
            !
            call self%perform_complex_synthesis(solution)

            !
            !==> Garbage collection
            !
            if ( xlmbda == ZERO ) then
                ! Nullify pointer
                nullify( iptr )
            else
                ! Release memory
                deallocate( iptr )
            end if
        end associate

    end subroutine invert_helmholtz



    subroutine get_gradient(self, scalar_function, &
        polar_gradient_component, azimuthal_gradient_component)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: scalar_function(:,:)
        real(wp),      intent(out)    :: polar_gradient_component(:,:)
        real(wp),      intent(out)    :: azimuthal_gradient_component(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_gradient'
        end if

        ! compute the (real) spherical harmonic coefficients
        call self%perform_scalar_analysis(scalar_function)

        ! compute gradient
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            v => polar_gradient_component, &
            w => azimuthal_gradient_component, &
            sqnn => self%vorticity_and_divergence_coefficients, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            br => self%workspace%real_polar_harmonic_coefficients, &
            bi => self%workspace%imaginary_polar_harmonic_coefficients, &
            cr => self%workspace%real_azimuthal_harmonic_coefficients, &
            ci => self%workspace%imaginary_azimuthal_harmonic_coefficients &
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
            call self%perform_vector_synthesis(v, w)
        end associate

    end subroutine get_gradient


    subroutine invert_gradient_from_spherical_components(self, &
        polar_source, azimuthal_source, solution )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: polar_source(:,:)
        real(wp),      intent(in)     :: azimuthal_source(:,:)
        real(wp),      intent(out)    :: solution(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in invert_gradient_from_spherical_components'
        end if

        ! Set vector spherical harmonic coefficients
        associate( &
            v => polar_source, &
            w => azimuthal_source &
            )
            call self%vector_analysis_from_spherical_components(v, w)
        end associate

        ! Invert gradient
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            f => solution, &
            sqnn => self%vorticity_and_divergence_coefficients, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            br => self%workspace%real_polar_harmonic_coefficients, &
            bi => self%workspace%imaginary_polar_harmonic_coefficients &
            )
            !
            !==> Initialize (real) coefficients
            !
            a = 0.0_wp
            b = 0.0_wp

            !
            !==> compute m=0 coefficients
            !
            do n=2, nlat
                a(1, n) = br(1, n)/sqnn(n)
                b(1, n) = bi(1, n)/sqnn(n)
            end do
            !
            !==> set upper limit for vector m subscript
            !
            associate( mmax => min(nlat, (nlon+1)/2) )
                !
                !==> compute m > 0 coefficients
                !
                do m=2, mmax
                    do n=m, nlat
                        a(m, n) = br(m, n)/sqnn(n)
                        b(m, n) = bi(m, n)/sqnn(n)
                    end do
                end do
            end associate
            !
            !==> Perform scalar synthesis
            !
            call self%perform_scalar_synthesis(f)
        end associate

    end subroutine invert_gradient_from_spherical_components



    subroutine get_vorticity_from_spherical_components(self, &
        polar_component, azimuthal_component, vorticity)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: polar_component(:,:)
        real(wp),      intent(in)     :: azimuthal_component(:,:)
        real(wp),      intent(out)    :: vorticity(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_vorticity_from_spherical_components'
        end if

        associate( &
            v => polar_component, &
            w => azimuthal_component, &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            vt => vorticity, &
            sqnn => self%vorticity_and_divergence_coefficients, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            cr => self%workspace%real_azimuthal_harmonic_coefficients, &
            ci => self%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            !
            !==> Perform vector analysis
            !
            call self%vector_analysis_from_spherical_components(v, w)
            !
            !==> Initialize (real) coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
            !
            !==> set upper limit for vector m subscript
            !
            associate( mmax => min(nlat, (nlon+1)/2) )
                !
                !==> compute m > 0 coefficients
                !
                do m=1, mmax
                    do n=m, nlat
                        a(m, n) = sqnn(n) * cr(m, n)
                        b(m, n) = sqnn(n) * ci(m, n)
                    end do
                end do
            end associate
            !
            !==> Compute vector harmonic synthesis
            !
            call self%perform_scalar_synthesis(vt)
        end associate


    end subroutine get_vorticity_from_spherical_components



    subroutine get_vorticity_from_vector_field(self, vector_field, vorticity)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: vector_field(:,:,:)
        real(wp),      intent(out)    :: vorticity(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        real(wp), allocatable :: polar_component(:,:)
        real(wp), allocatable :: azimuthal_component(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //'in get_vorticity_from_vector_field'
        end if

        !
        !==> Allocate memory
        !
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES &
            )
            allocate( polar_component(nlat, nlon) )
            allocate( azimuthal_component(nlat, nlon) )
        end associate

        !
        !==> Compute vorticity
        !
        associate( &
            F => vector_field, &
            v => polar_component, &
            w => azimuthal_component, &
            vort => vorticity &
            )
            !
            !==> Get spherical components
            !
            call self%unit_vectors%get_spherical_angle_components(F, v, w)
            !
            !==> Get vorticity from spherical components
            !
            call self%get_vorticity_from_spherical_components(v, w, vort)
        end associate

        !
        !==> Release memory
        !
        deallocate( polar_component )
        deallocate( azimuthal_component )

    end subroutine get_vorticity_from_vector_field



    subroutine invert_vorticity(self, source, polar_solution, azimuthal_solution)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: source(:,:)
        real(wp),      intent(out)    :: polar_solution(:,:)
        real(wp),      intent(out)    :: azimuthal_solution(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m !! Counters
        integer(ip) :: ityp_temp_save
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in invert_vorticity'
        end if

        !
        !==> compute the (real) spherical harmonic coefficients
        !
        call self%perform_scalar_analysis(source)

        !
        !==> Invert vorticity
        !
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            v => polar_solution, &
            w => azimuthal_solution, &
            sqnn => self%vorticity_and_divergence_coefficients, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            cr => self%workspace%real_azimuthal_harmonic_coefficients, &
            ci => self%workspace%imaginary_azimuthal_harmonic_coefficients, &
            isym => self%SCALAR_SYMMETRIES, &
            ityp => self%VECTOR_SYMMETRIES &
            )
            associate( mmax => min(nlat, (nlon+1)/2) )
                !
                !==> Initialize coefficients
                !
                cr = 0.0_wp
                ci = 0.0_wp

                !
                !==>  Compute m = 0 coefficients
                !
                do n = 2, nlat
                    cr(1, n) = a(1, n)/sqnn(n)
                    ci(1, n) = b(1, n)/sqnn(n)
                end do

                !
                !==> Compute m > 0 coefficients
                !
                do m=2, mmax
                    do n=m, nlat
                        cr(m, n) = a(m, n)/sqnn(n)
                        ci(m, n) = b(m, n)/sqnn(n)
                    end do
                end do

                !
                !==> Save old vector symmetry
                !
                ityp_temp_save = ityp

                !
                !==> Set symmetries for synthesis
                !
                select case (isym)
                    case (0)
                        ityp = 2
                    case (1)
                        ityp = 5
                    case(2)
                        ityp = 8
                end select

                !
                !==> Compute scalar synthesis
                !
                call self%perform_vector_synthesis(v, w)

                !
                !==> Reset old vector symmetry
                !
                ityp = ityp_temp_save
            end associate
        end associate

    end subroutine invert_vorticity



    subroutine get_divergence_from_vector_field(self, vector_field, divergence)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: vector_field (:,:,:)
        real(wp),      intent(out)    :: divergence (:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        real(wp), allocatable :: polar_component(:,:)
        real(wp), allocatable :: azimuthal_component(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //'in get_divergence_from_vector_field'
        end if

        !
        !==> Allocate memory
        !
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES &
            )
            allocate( polar_component(nlat, nlon) )
            allocate( azimuthal_component(nlat, nlon) )
        end associate

        !
        !==> Compute vorticity
        !
        associate( &
            F => vector_field, &
            v => polar_component, &
            w => azimuthal_component, &
            dv => divergence &
            )
            !
            !==> Get spherical components
            !
            call self%unit_vectors%get_spherical_angle_components(F, v, w)
            !
            !==> Get vorticity from spherical components
            !
            call self%get_divergence_from_spherical_components(v, w, dv)
        end associate

        !
        !==> Release memory
        !
        deallocate( polar_component )
        deallocate( azimuthal_component )

    end subroutine get_divergence_from_vector_field



    subroutine get_divergence_from_spherical_components(self, &
        polar_component, azimuthal_component, divergence )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: polar_component(:,:)
        real(wp),      intent(in)     :: azimuthal_component(:,:)
        real(wp),      intent(out)    :: divergence (:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_divergence_from_spherical_components'
        end if

        associate( &
            v => polar_component, &
            w => azimuthal_component, &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            dv => divergence, &
            sqnn => self%vorticity_and_divergence_coefficients, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            br => self%workspace%real_polar_harmonic_coefficients, &
            bi => self%workspace%imaginary_polar_harmonic_coefficients &
            )
            !
            !==> Perform vector analysis
            !
            call self%vector_analysis_from_spherical_components(v, w)
            !
            !==> Initialize (real) coefficients
            !
            a = 0.0_wp
            b = 0.0_wp
            !
            !==> set upper limit for vector m subscript
            !
            associate( mmax => min(nlat, (nlon+1)/2) )
                !
                !==> compute m > 0 coefficients
                !
                do m=1, mmax
                    do n=m, nlat
                        a(m, n) = -sqnn(n) * br(m, n)
                        b(m, n) = -sqnn(n) * bi(m, n)
                    end do
                end do
            end associate
            !
            !==> Compute vector harmonic synthesis
            !
            call self%perform_scalar_synthesis(dv)
        end associate
    end subroutine get_divergence_from_spherical_components



    subroutine invert_divergence(self, source, polar_solution, azimuthal_solution )
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: source(:,:)
        real(wp),      intent(out)    :: polar_solution(:,:)
        real(wp),      intent(out)    :: azimuthal_solution(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: n, m !! Counters
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in invert_divergence'
        end if

        ! Calculate the (real) scalar harmonic coefficients
        associate( div => source )

            call self%perform_scalar_analysis( div )

        end associate

        ! Invert gradient
        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            v => polar_solution, &
            w => azimuthal_solution, &
            sqnn => self%vorticity_and_divergence_coefficients, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            br => self%workspace%real_polar_harmonic_coefficients, &
            bi => self%workspace%imaginary_polar_harmonic_coefficients, &
            cr => self%workspace%real_azimuthal_harmonic_coefficients, &
            ci => self%workspace%imaginary_azimuthal_harmonic_coefficients &
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
            call self%perform_vector_synthesis(v, w)
        end associate

    end subroutine invert_divergence


    subroutine get_vorticity_and_divergence_from_velocities(self, &
        polar_component, azimuthal_component, vorticity, divergence)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere),         intent(inout)  :: self
        real(wp),              intent(in)     :: polar_component(:,:)
        real(wp),              intent(in)     :: azimuthal_component(:,:)
        real(wp),              intent(out)    :: vorticity(:,:)
        real(wp),              intent(out)    :: divergence(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_vorticity_and_divergence_from_velocities'
        end if

        associate( &
            v => polar_component, &
            w => azimuthal_component, &
            vt => vorticity, &
            dv => divergence &
            )
            !
            !==> Compute vorticity
            !
            call self%get_vorticity(v, w, vt)
            !
            !==> Compute divergence
            !
            call self%get_divergence(v, w, dv)
        end associate


    end subroutine get_vorticity_and_divergence_from_velocities


    subroutine get_velocities_from_vorticity_and_divergence_coefficients(self, &
        vort_spec, div_spec, polar_component, azimuthal_component)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        complex(wp),   intent(in)     :: vort_spec(:)
        complex(wp),   intent(in)     :: div_spec(:)
        real(wp),      intent(out)    :: polar_component(:,:)
        real(wp),      intent(out)    :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)           :: nm, n, m, i !! Counters
        real(wp), allocatable :: isqnn(:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_velocities_from_vorticity_and_divergence_coefficients'
        end if

        ! Associate various quantities
        associate( &
            v => polar_component, &
            w => azimuthal_component, &
            nlat => self%NUMBER_OF_LATITUDES, &
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            sqnn => self%vorticity_and_divergence_coefficients, &
            a => self%workspace%real_harmonic_coefficients, &
            b => self%workspace%imaginary_harmonic_coefficients, &
            br => self%workspace%real_polar_harmonic_coefficients, &
            bi => self%workspace%imaginary_polar_harmonic_coefficients, &
            cr => self%workspace%real_azimuthal_harmonic_coefficients, &
            ci => self%workspace%imaginary_azimuthal_harmonic_coefficients &
            )
            !
            !==> Allocate memory
            !
            allocate( isqnn(size(sqnn)) )

            isqnn(1) = 0.0_wp
            do n = 2, nlat
                isqnn(n) = 1.0_wp/sqnn(n)
            end do

            !
            !==> Initialize (real) scalar coefficients
            !
            a = 0.0_wp
            b = 0.0_wp

            !
            !==> Set (real) auxiliary coefficients
            !    from complex divergence harmonic coefficients
            !
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = -2.0_wp * real(div_spec(nm))
                    b(m, n) = -2.0_wp * aimag(div_spec(nm))
                end do
            end do

            !
            !==> Initialize polar vector coefficients
            !
            br = 0.0_wp
            bi = 0.0_wp

            !
            !==> Set polar vector coefficients
            !
            do n=1, nlat
                br(:, n) = isqnn(n)*a(:, n)
                bi(:, n) = isqnn(n)*b(:, n)
            end do

            !
            !==> Re-initialize (real) scalar coefficients
            !
            a = 0.0_wp
            b = 0.0_wp

            !
            !==> Set (real) auxiliary coefficients
            !    from complex vorticity harmonic coefficients
            !
            do m=1, ntrunc+1
                do n=m, ntrunc+1
                    nm = sum([(i, i=ntrunc+1, ntrunc-m+3, -1)])+n-m+1
                    a(m, n) = 2.0_wp * real(vort_spec(nm))
                    b(m, n) = 2.0_wp * aimag(vort_spec(nm))
                end do
            end do

            !
            !==> Initialize azimuthal vector coefficients
            !
            cr = 0.0_wp
            ci = 0.0_wp

            !
            !==> Set azimuthal vector coefficients
            !
            do n=1, nlat
                cr(:, n) = isqnn(n)*a(:, n)
                ci(:, n) = isqnn(n)*b(:, n)
            end do

            !
            !==> Compute vector harmonic sysnthesis to get components
            !
            call self%perform_vector_synthesis(v, w)
        end associate

        !
        !==> Release memory
        !
        deallocate( isqnn )

    end subroutine get_velocities_from_vorticity_and_divergence_coefficients

    subroutine get_velocities_from_vorticity_and_divergence(self, &
        vorticity, divergence, polar_component, azimuthal_component)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: vorticity(:,:)
        real(wp),      intent(in)     :: divergence(:,:)
        real(wp),      intent(out)    :: polar_component(:,:)
        real(wp),      intent(out)    :: azimuthal_component(:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        complex(wp), allocatable  :: vorticity_coefficients(:)
        complex(wp), allocatable  :: divergence_coefficients(:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_velocities_from_vorticity_and_divergence'
        end if

        associate( nm_dim => size(self%complex_spectral_coefficients) )
            !
            !==> Allocate memory
            !
            allocate( vorticity_coefficients(nm_dim) )
            allocate( divergence_coefficients(nm_dim) )

        end associate


        associate( &
            v => polar_component, &
            w => azimuthal_component, &
            vt_spec => vorticity_coefficients, &
            dv_spec => divergence_coefficients, &
            vt => vorticity, &
            dv => divergence &
            )
            !
            !==> Compute complex spectral coefficients
            !
            call self%analyze_into_complex_spectral_coefficients(vt, vt_spec)
            call self%analyze_into_complex_spectral_coefficients(dv, dv_spec)
            !
            !==> Compute velocities
            !
            call self%get_velocities_from_vorticity_and_divergence_coefficients( &
                vt_spec, dv_spec, v, w)
        end associate


    end subroutine get_velocities_from_vorticity_and_divergence

    subroutine compute_angular_momentum(self, scalar_function, angular_momentum)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        real(wp),      intent(in)     :: scalar_function(:,:)
        real(wp),      intent(out)    :: angular_momentum(:,:,:)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)           :: k, l   !! Counters
        real(wp), allocatable :: polar_gradient_component(:,:)
        real(wp), allocatable :: azimuthal_gradient_component(:,:)
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in compute_angular_momentum'
        end if

        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES &
            )
            !
            !==>  Allocate memory
            !
            allocate( polar_gradient_component(nlat, nlon) )
            allocate( azimuthal_gradient_component(nlat, nlon) )
        end associate


        associate( &
            f => scalar_function, &
            grad_theta => polar_gradient_component, &
            grad_phi => azimuthal_gradient_component &
            )
            !
            !==> Calculate the spherical surface gradient components
            !
            call self%get_gradient(f, grad_theta, grad_phi)

        end associate


        associate( &
            nlat => self%NUMBER_OF_LATITUDES, &
            nlon => self%NUMBER_OF_LONGITUDES, &
            R => angular_momentum &
            )
            do l = 1, nlon
                do k = 1, nlat
                    associate( &
                        theta => self%unit_vectors%polar(k, l), &
                        phi => self%unit_vectors%azimuthal(k, l), &
                        grad_theta => polar_gradient_component(k, l), &
                        grad_phi => azimuthal_gradient_component(k, l) &
                        )
                        !
                        !==> Calculate the rotation operator applied to a scalar function
                        !
                        R(:, k, l) = phi * grad_theta - theta * grad_phi

                    end associate
                end do
            end do
        end associate

        !
        !==> Release memory
        !
        deallocate( polar_gradient_component )
        deallocate( azimuthal_gradient_component )

    end subroutine compute_angular_momentum
    


    function get_index(self, n, m) result (return_value)
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
        ! integer(ip), dimension ((MTRUNC+1)*(MTRUNC+2)/2) :: indxm, indxn
        ! indxm = [((m, n=m, MTRUNC), m=0, MTRUNC)]
        ! indxn = [((n, n=m, MTRUNC), m=0, MTRUNC)]
        !
        ! Conversely, the index nm as a function of m and n is:
        ! nm = sum([(i, i=MTRUNC+1, MTRUNC-m+2, -1)])+n-m+1
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        integer(ip),   intent(in)     :: n
        integer(ip),   intent(in)     :: m
        integer(ip)                    :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: i !! Counter
        !----------------------------------------------------------------------

        ! Initialize return value
        return_value = -1

        associate( ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT )

            if ( m <= n .and. max(n, m) <= ntrunc ) then
                return_value = sum ([(i, i=ntrunc+1, ntrunc-m+2, -1)]) + n-m+1
            end if

        end associate

    end function get_index


    function get_coefficient(self, n, m) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        integer(ip),   intent(in)     :: n
        integer(ip),   intent(in)     :: m
        complex(wp)                    :: return_value
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (.not.self%initialized) then
            error stop 'Uninitialized object of class(Sphere): '&
                //' in get_coefficient'
        end if

        associate( &
            ntrunc => self%TRIANGULAR_TRUNCATION_LIMIT, &
            nm  => self%get_index(n, m), &
            nm_conjg => self%get_index(n, -m), &
            psi => self%complex_spectral_coefficients &
            )

            if (m < 0 .and. nm_conjg > 0) then
                return_value = ( (-1.0_wp)**(-m) ) * conjg(psi(nm_conjg))
            else if (nm > 0) then
                return_value = psi(nm)
            else
                return_value = 0.0_wp
            end if

        end associate

    end function get_coefficient


    subroutine get_scalar_symmetries(self, isym)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        integer(ip),   intent(in)     :: isym
        !----------------------------------------------------------------------

        select case (isym)
            case (0:2)
                self%SCALAR_SYMMETRIES = isym
            case default
                error stop 'Object of class(Sphere) in get_scalar_symmetries: '&
                    //'Optional argument isym must be either '&
                    //'0, 1, or 2 (default isym = 0)'
        end select

    end subroutine get_scalar_symmetries
    
    subroutine get_vector_symmetries(self, itype)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class(Sphere), intent(inout)  :: self
        integer(ip),   intent(in)     :: itype
        !----------------------------------------------------------------------

        select case (itype)
            case (0:8)
                self%VECTOR_SYMMETRIES = itype
            case default
                error stop 'Object of class(Sphere) in get_vector_symmetries: '&
                    //'Optional argument itype must be either '&
                    //'0, 1, 2, ..., 8 (default itype = 0)'
        end select

    end subroutine get_vector_symmetries

end module type_Sphere
