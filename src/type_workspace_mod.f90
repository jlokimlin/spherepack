!*****************************************************************************************
!
! Purpose:
!
! Defines the derived data type "workspace_t" required to
! invoke the spherepack library
!
!*****************************************************************************************
!
module type_workspace_mod

    use, intrinsic :: iso_fortran_env, only: &
        WP     => REAL64, &
        IP     => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: workspace_t

    !--------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------- ---------------------------------------
    character (len=200) :: error_message            !! Probably long enough
    integer (IP)        :: allocate_status          !! To check allocation status
    integer (IP)        :: deallocate_status        !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, public :: workspace_t

        ! All components are public unless stated otherwise

        !---------------------------------------------------------------------------------
        ! Initialization flag
        !---------------------------------------------------------------------------------
        logical                         :: initialized = .false.
        !---------------------------------------------------------------------------------
        ! Workspace arrays for Legendre transform
        !---------------------------------------------------------------------------------
        real (WP), allocatable          :: work(:)
        real (WP), allocatable, private :: dwork(:)
        !---------------------------------------------------------------------------------
        ! Workspace arrays for scalar transform - GAU
        !---------------------------------------------------------------------------------
        real (WP), allocatable          :: wshags(:)
        real (WP), allocatable          :: wshsgs(:)
        !---------------------------------------------------------------------------------
        ! Workspace arrays for scalar transform - REG
        !---------------------------------------------------------------------------------
        real (WP), allocatable          :: wshaes(:)
        real (WP), allocatable          :: wshses(:)
        !---------------------------------------------------------------------------------
        ! Workspace arrays for vector transform - GAU
        !---------------------------------------------------------------------------------
        real (WP), allocatable          :: wvhags(:)
        real (WP), allocatable          :: wvhsgs(:)
        !---------------------------------------------------------------------------------
        ! Workspace arrays for vector transform - REG
        !---------------------------------------------------------------------------------
        real (WP), allocatable         :: wvhaes(:)
        real (WP), allocatable         :: wvhses(:)
        !---------------------------------------------------------------------------------
        ! Scalar transform coefficients
        !---------------------------------------------------------------------------------
        real (WP), pointer            :: real_harmonic_coefficients(:,:)      => null()
        real (WP), pointer            :: imaginary_harmonic_coefficients(:,:) => null()
        !---------------------------------------------------------------------------------
        ! Vector transform coefficients
        !---------------------------------------------------------------------------------
        real (WP), pointer            :: real_polar_harmonic_coefficients(:,:)          => null()
        real (WP), pointer            :: imaginary_polar_harmonic_coefficients(:,:)     => null()
        real (WP), pointer            :: real_azimuthal_harmonic_coefficients(:,:)      => null()
        real (WP), pointer            :: imaginary_azimuthal_harmonic_coefficients(:,:) => null()
        !---------------------------------------------------------------------------------

    contains

        ! All methods are private unless stated otherwise
        private

        !---------------------------------------------------------------------------------
        procedure                          :: Assert_initialized
        procedure, nopass                  :: Get_lwork
        procedure, nopass                  :: Get_ldwork
        !---------------------------------------------------------------------------------
        ! Public methods
        !---------------------------------------------------------------------------------
        procedure, non_overridable, public :: Create
        procedure, non_overridable, public :: Destroy
        !---------------------------------------------------------------------------------
        ! Private methods for gaussian grids
        !---------------------------------------------------------------------------------
        procedure, nopass                  :: Get_lshags
        procedure, nopass                  :: Get_lshsgs
        procedure, nopass                  :: Get_lvhags
        procedure, nopass                  :: Get_lvhsgs
        procedure                          :: Initialize_scalar_analysis_gau
        procedure                          :: Initialize_scalar_synthesis_gau
        procedure                          :: Initialize_scalar_transform_gau
        procedure                          :: Initialize_vector_analysis_gau
        procedure                          :: Initialize_vector_synthesis_gau
        procedure                          :: Initialize_vector_transform_gau
        !---------------------------------------------------------------------------------
        ! Private methods for regular (equally-spaced) grids
        !---------------------------------------------------------------------------------
        procedure, nopass                  :: Get_lshaes
        procedure, nopass                  :: Get_lshses
        procedure, nopass                  :: Get_lvhaes
        procedure, nopass                  :: Get_lvhses
        procedure                          :: Initialize_scalar_analysis_reg
        procedure                          :: Initialize_scalar_synthesis_reg
        procedure                          :: Initialize_scalar_transform_reg
        procedure                          :: Initialize_vector_analysis_reg
        procedure                          :: Initialize_vector_synthesis_reg
        procedure                          :: Initialize_vector_transform_reg
        !---------------------------------------------------------------------------------
        ! Finalizer
        !---------------------------------------------------------------------------------
        final                              :: Finalize_workspace
        !---------------------------------------------------------------------------------

    end type workspace_t

contains
    !
    !*****************************************************************************************
    !
    subroutine Create( this, nlat, nlon, grid_type )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out)       :: this
        integer (IP),        intent (in)           :: nlat
        integer (IP),        intent (in)           :: nlon
        character (len=*),   intent (in), optional :: grid_type
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization status
        !--------------------------------------------------------------------------------

        if ( this%initialized ) then

            write( stderr, '(A)' ) 'ERROR: TYPE(workspace_t)'
            write( stderr, '(A)')  ' You must destroy object before re-instantiating'

        end if

        !--------------------------------------------------------------------------------
        ! Allocate desired workspace arrays
        !--------------------------------------------------------------------------------

        ! Check if optional grid type argument is present
        if ( present( grid_type ) ) then

            ! Check if the grid type is equally spaced
            if ( grid_type .eq. 'REG' ) then

                ! Set up transforms for regular grids
                call this%Initialize_scalar_transform_reg( nlat, nlon )
                call this%Initialize_vector_transform_reg( nlat, nlon )

            ! Check if the grid type is gaussian
            else if ( grid_type .eq. 'GAU' ) then

                ! Set up transforms for gaussian grids
                call this%Initialize_scalar_transform_gau( nlat, nlon )
                call this%Initialize_vector_transform_gau( nlat, nlon )

            else

                ! Handle invalid grid type
                write( stderr, '(A)' ) 'ERROR: TYPE(workspace_t)'
                write( stderr, '(A)' ) 'ERROR: optional argument grid_type = ', grid_type
                write( stderr, '(A)' ) 'must be either REG or GAU (default GAU)'
                stop

            end if
        else

            ! Set up default transforms for gaussian grids
            call this%Initialize_scalar_transform_gau( nlat, nlon )
            call this%Initialize_vector_transform_gau( nlat, nlon )

        end if

        !--------------------------------------------------------------------------------
        ! Set initialization status
        !--------------------------------------------------------------------------------

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
        class (workspace_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization status
        !--------------------------------------------------------------------------------

        if ( .not. this%initialized ) return

        !--------------------------------------------------------------------------------
        ! Deallocate workspace arrays for Legendre transform
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%work ) ) then

            ! Deallocate array
            deallocate ( &
                this%work, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WORK failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Deallocate workspace arrays for scalar transform - GAU
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%wshags ) ) then

            ! Deallocate array
            deallocate ( &
                this%wshags, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHAGS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wshsgs ) ) then
            deallocate ( &
                this%wshsgs, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHSGS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Deallocate workspace arrays for scalar transform - REG
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%wshaes ) ) then

            ! Deallocate array
            deallocate ( &
                this%wshaes, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHAES failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wshses ) ) then

            ! Deallocate array
            deallocate ( &
                this%wshses, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHSES failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Deallocate workspace arrays for vector transform - GAU
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%wvhags ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhags, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHAGS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wvhsgs ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhsgs, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHSGS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Deallocate workspace arrays for vector transform - REG
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%wvhaes ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhaes, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHAES failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wvhses ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhses, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHSES failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Clean up scalar transform pointers
        !--------------------------------------------------------------------------------

        ! Check if pointer is associated
        if ( associated(this%real_harmonic_coefficients) ) then

            ! Deallocate pointer
            deallocate( &
                this%real_harmonic_coefficients, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating REAL_HARMONIC_COEFFICIENTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

            ! Nullify pointer
            nullify( this%real_harmonic_coefficients )

        end if

        ! Check if pointer is associated
        if ( associated(this%imaginary_harmonic_coefficients) ) then

            ! Deallocate pointer
            deallocate( &
                this%imaginary_harmonic_coefficients, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating IMAGINARY_HARMONIC_COEFFICIENTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

            ! Nullify pointer
            nullify( this%imaginary_harmonic_coefficients )

        end if

        !--------------------------------------------------------------------------------
        ! Clean up vector transform pointers
        !--------------------------------------------------------------------------------

        ! Check if pointer is associated
        if ( associated(this%real_polar_harmonic_coefficients) ) then

            ! Deallocate pointer
            deallocate( &
                this%real_polar_harmonic_coefficients, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating REAL_POLAR_HARMONIC_COEFFICIENTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

            ! Nullify pointer
            nullify( this%real_polar_harmonic_coefficients )

        end if

        ! Check if pointer is associated
        if ( associated(this%imaginary_polar_harmonic_coefficients) ) then

            ! Deallocate pointer
            deallocate( &
                this%imaginary_polar_harmonic_coefficients, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating IMAGINARY_POLAR_HARMONIC_COEFFICIENTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

            ! Nullify pointer
            nullify( this%imaginary_polar_harmonic_coefficients )

        end if

        ! Check if pointer is associated
        if ( associated(this%real_azimuthal_harmonic_coefficients) ) then

            ! Deallocate pointer
            deallocate( &
                this%real_azimuthal_harmonic_coefficients, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating REAL_AZIMUTHAL_HARMONIC_COEFFICIENTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

            ! Nullify pointer
            nullify( this%real_azimuthal_harmonic_coefficients )

        end if

        ! Check if pointer is associated
        if ( associated(this%imaginary_azimuthal_harmonic_coefficients) ) then

            ! Deallocate pointer
            deallocate( &
                this%imaginary_azimuthal_harmonic_coefficients, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating IMAGINARY_AZIMUTHAL_HARMONIC_COEFFICIENTS failed in DESTROY'
                write( stderr, '(A)' ) trim( error_message )

            end if

            ! Nullify pointer
            nullify( this%imaginary_azimuthal_harmonic_coefficients )

        end if

        !--------------------------------------------------------------------------------
        ! Reset initialization status
        !--------------------------------------------------------------------------------

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
        class (workspace_t), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        ! Check status
        if ( this%initialized ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'You must re-instantiate object with a call to DESTROY before calling methods'

        end if

    end subroutine Assert_initialized
    !
    !*****************************************************************************************
    !
    ! Private methods for gaussian grids
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_analysis_gau( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        !
        !--------------------------------------------------------------------------------
        !
        !    Documentation: SPHEREPACK 3.2
        !
        !     subroutine shagsi(nlat,nlon,wshags,lshags,work,lwork,dwork,ldwork,
        !                      ierror)
        !
        !     subroutine shagsi initializes the array wshags which can then
        !     be used repeatedly by subroutines shags. it precomputes
        !     and stores in wshags quantities such as gaussian weights,
        !     legendre polynomial coefficients, and fft trigonometric tables.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0,pi) and are compu
        !            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
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
        !     wshags an array which must be initialized by subroutine shagsi.
        !            once initialized, wshags can be used repeatedly by shags
        !            as long as nlat and nlon remain unchanged.  wshags must
        !            not be altered between calls of shags.
        !
        !     lshags the dimension of the array wshags as it appears in the
        !            program that calls shags. define
        !
        !               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
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
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls shagsi. lwork must be at least
        !            4*nlat*(nlat+2)+2 in the routine calling shagsi
        !
        !     dwork   a double precision work array that does not have to be saved.
        !
        !     ldwork  the length of dwork in the calling routine.  ldwork must
        !             be at least nlat*(nlat+4)
        !
        !     output parameter
        !
        !     wshags an array which must be initialized before calling shags or
        !            once initialized, wshags can be used repeatedly by shags or
        !            as long as nlat and nlon remain unchanged.  wshags must not
        !            altered between calls of shasc.
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lshags
        !            = 4  error in the specification of lwork
        !            = 5  error in the specification of ldwork
        !            = 6  failure in gaqd to compute gaussian points
        !                 (due to failure in eigenvalue routine)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: error_flag
        integer (IP) :: lwork, ldwork, lshags
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Compute dimensions of various workspace arrays
        !--------------------------------------------------------------------------------

        lwork  = Get_lwork( nlat, nlon )
        ldwork = Get_ldwork( nlat )
        lshags = Get_lshags( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%work ) ) then

            ! Deallocate array
            deallocate ( &
                this%work, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WORK failed in INITIALIZE_SCALAR_ANALYSIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_SCALAR_ANALYSIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wshags ) ) then

            ! Deallocate array
            deallocate ( &
                this%wshags, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHAGS failed in INITIALIZE_SCALAR_ANALYSIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        allocate ( &
            this%work(   1:lwork ), &
            this%dwork(  1:ldwork ), &
            this%wshags( 1:lshags ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_SCALAR_ANALYSIS_GAU'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2 routine
        !--------------------------------------------------------------------------------

        associate( &
            wshags => this%wshags, &
            work   => this%work, &
            dwork  => this%dwork,&
            ierror => error_flag &
            )

            call Shagsi( nlat, nlon, wshags, lshags, work, lwork, dwork, ldwork, &
                ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_SCALAR_ANALYSIS_GAU'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for SHAGS'
                write( stderr, '(A)') 'size( this%SHAGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for WORK'
                write( stderr, '(A)') 'size( this%WORK )'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Error in the specification of extent for DWORK'
                write( stderr, '(A)') 'size( this%DWORK )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Error in call to GAQD to compute gaussian points'
                write( stderr, '(A)') '(due to failure in eigenvalue routine)'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_scalar_analysis_gau
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_synthesis_gau( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic synthesis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        !
        !--------------------------------------------------------------------------------
        !
        !     Documentation: SPHEREPACK 3.2
        !
        !     subroutine shsgsi(nlat,nlon,wshsgs,lshsgs,work,lwork,dwork,ldwork,
        !                     ierror)
        !
        !     subroutine shsgsi initializes the array wshsgs which can then
        !     be used repeatedly by subroutines shsgs. it precomputes
        !     and stores in wshsgs quantities such as gaussian weights,
        !     legendre polynomial coefficients, and fft trigonometric tables.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0,pi) and are compu
        !            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
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
        !     wshsgs an array which must be initialized by subroutine shsgsi.
        !            once initialized, wshsgs can be used repeatedly by shsgs
        !            as long as nlat and nlon remain unchanged.  wshsgs must
        !            not be altered between calls of shsgs.
        !
        !     lshsgs the dimension of the array wshsgs as it appears in the
        !            program that calls shsgs. define
        !
        !               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
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
        !     work   a real work space which need not be saved
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls shsgsi. lwork must be at least
        !            4*nlat*(nlat+2)+2 in the routine calling shsgsi
        !
        !     dwork   a double precision work array that does not have to be saved.
        !
        !     ldwork  the length of dwork in the calling routine.  ldwork must
        !             be at least nlat*(nlat+4)
        !
        !     output parameter
        !
        !     wshsgs an array which must be initialized before calling shsgs or
        !            once initialized, wshsgs can be used repeatedly by shsgs or
        !            as long as nlat and nlon remain unchanged.  wshsgs must not
        !            altered between calls of shsgs.
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lshsgs
        !            = 4  error in the specification of lwork
        !            = 5  error in the specification of ldwork
        !            = 6  failure in gaqd to compute gaussian points
        !                 (due to failure in eigenvalue routine)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: error_flag
        integer (IP) :: lwork, ldwork, lshsgs
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Compute dimensions of various workspace arrays
        !--------------------------------------------------------------------------------

        lwork  = Get_lwork( nlat, nlon )
        ldwork = Get_ldwork( nlat )
        lshsgs = Get_lshsgs( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%work ) ) then

            ! Deallocate array
            deallocate ( &
                this%work, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WORK failed in INITIALIZE_SCALAR_SYNTHESIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_SCALAR_SYNTHESIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wshsgs ) ) then
            deallocate ( &
                this%wshsgs, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHSGS failed in INITIALIZE_SCALAR_SYNTHESIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        ! Allocate arrays
        allocate ( &
            this%work(   1:lwork ), &
            this%dwork(  1:ldwork ), &
            this%wshsgs( 1:lshsgs ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_SCALAR_SYNTHESIS_GAU'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            wshsgs => this%wshsgs, &
            work   => this%work, &
            dwork  => this%dwork, &
            ierror => error_flag &
            )

            call Shsgsi( nlat, nlon, wshsgs, lshsgs, work, lwork, dwork, ldwork, &
                ierror)

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_SCALAR_SYNTHESIS_GAU'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for SHSGS'
                write( stderr, '(A)') 'size( this%SHSGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for WORK'
                write( stderr, '(A)') 'size( this%WORK )'

            else if ( error_flag == 5 ) then

                write( stderr, '(A)') 'Error in the specification of extent for DWORK'
                write( stderr, '(A)') 'size( this%DWORK )'

            else if ( error_flag == 6 ) then

                write( stderr, '(A)') 'Error in call to GAQD to compute gaussian points'
                write( stderr, '(A)') '(due to failure in eigenvalue routine)'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_scalar_synthesis_gau
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_transform_gau( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays and pointers to perform the
        !  (real) scalar harmonic transform on a gaussian grid
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Set up scalar analysis
        !--------------------------------------------------------------------------------

        call this%Initialize_scalar_analysis_gau( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Set up scalar synthesis
        !--------------------------------------------------------------------------------

        call this%Initialize_scalar_synthesis_gau( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Allocate pointers for the (real) scalar transform
        !--------------------------------------------------------------------------------

        allocate ( &
            this%real_harmonic_coefficients( 1:nlat, 1:nlat ), &
            this%imaginary_harmonic_coefficients( 1:nlat, 1:nlat ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Pointer allocation failed in INITIALIZE_SCALAR_TRANSFORM_GAU'
            write( stderr, '(A)' )  trim( error_message )

        end if

    end subroutine Initialize_scalar_transform_gau
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_analysis_gau( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#vhags.html
        !
        !--------------------------------------------------------------------------------
        !
        ! Documentation: SPHEREPACK 3.2
        !
        !     subroutine vhagsi(nlat,nlon,wvhags,lvhags,work,lwork,ierror)
        !
        !     subroutine vhagsi initializes the array wvhags which can then be
        !     used repeatedly by subroutine vhags until nlat or nlon is changed.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0,pi) and are computed
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
        !     lvhags the dimension of the array wvhags as it appears in the
        !            program that calls vhagsi.  define
        !
        !               l1 = min0(nlat,nlon/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lvhags must be at least
        !
        !               3*nlat*(nlat+1)+2  (required by vhagsi)
        !
        !     dwork  a double precision work space that does not need to be saved
        !
        !     ldwork the dimension of the array dwork as it appears in the
        !            program that calls vhagsi. ldwork must be at least
        !
        !                   (3*nlat*(nlat+3)+2)/2
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     wvhags an array which is initialized for use by subroutine vhags.
        !            once initialized, wvhags can be used repeatedly by vhags
        !            as long as nlat and nlon remain unchanged.  wvhags must not
        !            be altered between calls of vhags.
        !
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lvhags
        !            = 4  error in the specification of ldwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP), intent (in)            :: nlat
        integer (IP), intent (in)            :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: error_flag
        integer (IP)    :: ldwork, lvhags
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Compute dimensions of various workspace arrays
        !--------------------------------------------------------------------------------

        ldwork = Get_ldwork( nlat )
        lvhags = Get_lvhags( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_SCALAR_ANALYSIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wvhags ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhags, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHAGS failed in INITIALIZE_SCALAR_ANALYSIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        allocate ( &
            this%dwork(  1:ldwork ), &
            this%wvhags( 1:lvhags ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_SCALAR_ANALYSIS_GAU'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2 routine
        !--------------------------------------------------------------------------------

        associate( &
            wvhags => this%wvhags, &
            dwork  => this%dwork,&
            ierror => error_flag &
            )

            call Vhagsi( nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_VECTOR_ANALYSIS_GAU'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for WVHAGS'
                write( stderr, '(A)') 'size( this%WVHAGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for DWORK'
                write( stderr, '(A)') 'size( this%DWORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_vector_analysis_gau
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_synthesis_gau( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#vhsgsi.html
        !
        !--------------------------------------------------------------------------------
        !
        ! Documentation: SPHEREPACK 3.2
        !
        !     subroutine vhsgsi(nlat,nlon,wvhsgs,lvhsgs,work,lwork,ierror)
        !
        !     subroutine vhsgsi initializes the array wvhsgs which can then be
        !     used repeatedly by subroutine vhsgs until nlat or nlon is changed.
        !
        !     input parameters
        !
        !     nlat   the number of points in the gaussian colatitude grid on the
        !            full sphere. these lie in the interval (0,pi) and are computed
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
        !     lvhsgs the dimension of the array wvhsgs as it appears in the
        !            program that calls vhsgs. define
        !
        !               l1 = min0(nlat,nlon/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
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
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls vhsgsi. lwork must be at least
        !
        !               4*nlat*nlon*(2*nt+1)
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     wvhsgs an array which is initialized for use by subroutine vhsgs.
        !            once initialized, wvhsgs can be used repeatedly by vhsgs
        !            as long as nlat and nlon remain unchanged.  wvhsgs must not
        !            be altered between calls of vhsgs.
        !
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lvhsgs
        !            = 4  error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: error_flag
        integer (IP)    :: ldwork, lvhsgs
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Compute dimensions of various workspace arrays
        !--------------------------------------------------------------------------------

        ldwork = Get_ldwork( nlat )
        lvhsgs = Get_lvhsgs( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_VECTOR_SYNTHESIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wvhsgs ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhsgs, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHSGS failed in INITIALIZE_VECTOR_SYNTHESIS_GAU'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        allocate ( &
            this%dwork(  1:ldwork ), &
            this%wvhsgs( 1:lvhsgs ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in in INITIALIZE_VECTOR_SYNTHESIS_GAU'
            write( stderr, '(A)' ) trim( error_message )

        end if


        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------


        associate( &
            wvhsgs => this%wvhsgs, &
            dwork  => this%dwork, &
            ierror => error_flag &
            )

            call Vhsgsi( nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_VECTOR_ANALYSIS_GAU'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for WVHSGS'
                write( stderr, '(A)') 'size( this%WVHSGS )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for DWORK'
                write( stderr, '(A)') 'size( this%DWORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_vector_synthesis_gau
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_transform_gau( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Sets the various workspace arrays and pointers
        !  required for the (real) vector harmonic transform
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Set up vector analysis
        !--------------------------------------------------------------------------------

        call this%Initialize_vector_analysis_gau( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Set up vector analysis
        !--------------------------------------------------------------------------------

        call this%Initialize_vector_synthesis_gau( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Allocate pointers for the (real) vector harmonic transform coefficients
        !--------------------------------------------------------------------------------

        allocate ( &
            this%real_polar_harmonic_coefficients( 1:nlat, 1:nlat ), &
            this%imaginary_polar_harmonic_coefficients( 1:nlat, 1:nlat ), &
            this%real_azimuthal_harmonic_coefficients( 1:nlat, 1:nlat ), &
            this%imaginary_azimuthal_harmonic_coefficients( 1:nlat, 1:nlat ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocating pointers failed in INITIALIZE_VECTOR_TRANSFORM_GAU'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end subroutine Initialize_vector_transform_gau
    !
    !*****************************************************************************************
    !
    ! Private methods for regular (equally_spaced) grids
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_analysis_reg( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a regular (equally-spaced) grid.
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shaesi.html
        !
        !--------------------------------------------------------------------------------
        !
        !  Documentation: SPHEREPACK 3.2
        !
        !     subroutine shaesi(nlat,nlon,wshaes,lshaes,work,lwork,dwork,
        !                      ldwork,ierror)
        !
        !     subroutine shaesi initializes the array wshaes which can then
        !     be used repeatedly by subroutine shaes
        !
        !     input parameters
        !
        !     nlat   the number of colatitudes on the full sphere including the
        !            poles. for example, nlat = 37 for a five degree grid.
        !            nlat determines the grid increment in colatitude as
        !            pi/(nlat-1).  if nlat is odd the equator is located at
        !            grid point i=(nlat+1)/2. if nlat is even the equator is
        !            located half way between points i=nlat/2 and i=nlat/2+1.
        !            nlat must be at least 3. note: on the half sphere, the
        !            number of grid points in the colatitudinal direction is
        !            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than or equal to 4. the efficiency of the computation is
        !            improved when nlon is a product of small prime numbers.
        !
        !     lshaes the dimension of the array wshaes as it appears in the
        !            program that calls shaesi. define
        !
        !               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lshaes must be at least
        !
        !               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
        !
        !     work   a real   work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls shaesi.  define
        !
        !               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lwork must be at least
        !
        !               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
        !
        !
        !     dwork  a double precision work array that does not have to be saved.
        !
        !     ldwork the dimension of the array dwork as it appears in the
        !            program that calls shaesi.  ldwork must be at least nlat+1
        !
        !
        !     output parameters
        !
        !     wshaes an array which is initialized for use by subroutine shaes.
        !            once initialized, wshaes can be used repeatedly by shaes
        !            as long as nlon and nlat remain unchanged.  wshaes must
        !            not be altered between calls of shaes.
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lshaes
        !            = 4  error in the specification of lwork
        !            = 5  error in the specification of ldwork
        !
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: error_flag
        integer (IP) :: lwork, ldwork, lshaes
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Compute dimensions of various workspace arrays
        !--------------------------------------------------------------------------------

        lwork  = this%Get_lwork(  nlat, nlon )
        ldwork = this%Get_ldwork( nlat )
        lshaes = this%Get_lshaes( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%work ) ) then

            ! Deallocate array
            deallocate ( &
                this%work, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WORK failed in INITIALIZE_SCALAR_ANALYSIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_SCALAR_ANALYSIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wshaes ) ) then

            ! Deallocate array
            deallocate ( &
                this%wshaes, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHAES failed in INITIALIZE_SCALAR_ANALYSIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        allocate ( &
            this%work(   1:lwork ), &
            this%dwork(  1:ldwork ), &
            this%wshaes( 1:lshaes ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_SCALAR_ANALYSIS_REG'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            wshaes => this%wshaes, &
            work   => this%work, &
            dwork  => this%dwork,  &
            ierror => error_flag &
            )

            call Shaesi( nlat, nlon, wshaes, lshaes, work, lwork, dwork, ldwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_SCALAR_ANALYSIS_REG'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for WSHAES'
                write( stderr, '(A)') 'size( this%WSHAES )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for WORK'
                write( stderr, '(A)') 'size( this%WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_scalar_analysis_reg
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_synthesis_reg( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a regular (equally-spaced) grid.
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shsesi.html
        !
        !--------------------------------------------------------------------------------
        !
        !  Documentation: SPHEREPACK 3.2
        !
        !     subroutine shsesi(nlat,nlon,wshses,lshses,work,lwork,ierror)
        !
        !     subroutine shsesi initializes the array wshses which can then
        !     be used repeatedly by subroutine shses.
        !
        !     input parameters
        !
        !     nlat   the number of colatitudes on the full sphere including the
        !            poles. for example, nlat = 37 for a five degree grid.
        !            nlat determines the grid increment in colatitude as
        !            pi/(nlat-1).  if nlat is odd the equator is located at
        !            grid point i=(nlat+1)/2. if nlat is even the equator is
        !            located half way between points i=nlat/2 and i=nlat/2+1.
        !            nlat must be at least 3. note: on the half sphere, the
        !            number of grid points in the colatitudinal direction is
        !            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than or equal to 4. the efficiency of the computation is
        !            improved when nlon is a product of small prime numbers.
        !
        !     lshses the dimension of the array wshses as it appears in the
        !            program that calls shsesi. define
        !
        !               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lshses must be at least
        !
        !               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls shsesi.  define
        !
        !               l1 = min0(nlat,(nlon+2)/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lwork must be at least
        !
        !               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
        !
        !     output parameters
        !
        !     wshses an array which is initialized for use by subroutine shses.
        !            once initialized, wshses can be used repeatedly by shses
        !            as long as nlon and nlat remain unchanged.  wshses must
        !            not be altered between calls of shses.
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lshses
        !            = 4  error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) :: error_flag
        integer (IP) :: lwork, ldwork, lshses
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Set up various workspace dimensions
        !--------------------------------------------------------------------------------

        lwork  = this%Get_lwork(  nlat, nlon )
        ldwork = this%Get_ldwork( nlat )
        lshses = this%Get_lshses( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%work ) ) then

            ! Deallocate array
            deallocate ( &
                this%work, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WORK failed in INITIALIZE_SCALAR_ANALYSIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_SCALAR_ANALYSIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wshses ) ) then

            ! Deallocate array
            deallocate ( &
                this%wshses, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WSHSES failed in INITIALIZE_SCALAR_SYNTHESIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        allocate ( &
            this%work(   1:lwork ), &
            this%dwork(  1:ldwork ), &
            this%wshses( 1:lshses ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_SCALAR_SYNTHESIS_REG'
            write( stderr, '(A)' ) trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            wshses => this%wshses, &
            work   => this%work, &
            dwork  => this%dwork, &
            ierror => error_flag &
            )

            call Shsesi( nlat, nlon, wshses, lshses, work, lwork, dwork, ldwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_SCALAR_SYNTHESIS_REG'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for WSHSES'
                write( stderr, '(A)') 'size( this%WSHSES )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for WORK'
                write( stderr, '(A)') 'size( this%WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_scalar_synthesis_reg
    !
    !*****************************************************************************************
    !
    subroutine Initialize_scalar_transform_reg( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays and pointers to perform the
        !  (real) scalar harmonic transform on a regular (equally-spaced) grid.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Set up scalar analysis - REG
        !--------------------------------------------------------------------------------

        call this%Initialize_scalar_analysis_reg( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Set up scalar synthesis - REG
        !--------------------------------------------------------------------------------

        call this%Initialize_scalar_synthesis_reg( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Allocate pointers for the (real) scalar harmonic transform - REG
        !--------------------------------------------------------------------------------

        allocate ( &
            this%real_harmonic_coefficients(      1:nlat, 1:nlat ), &
            this%imaginary_harmonic_coefficients( 1:nlat, 1:nlat ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Pointer allocation failed in INITIALIZE_SCALAR_TRANSFORM_REG'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end subroutine Initialize_scalar_transform_reg
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_analysis_reg( this, nlat, nlon )
        !
        !< Purpose:
        !
        ! Sets the various workspace arrays required for
        ! (real) vector harmonic analysis on a regular (equally-spaced) grid
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhaesi.html
        !
        !--------------------------------------------------------------------------------
        !
        !     Documentation: SPHEREPACK 3.2
        !
        !     subroutine vhaesi(nlat,nlon,wvhaes,lvhaes,work,lwork,ierror)
        !
        !     subroutine vhaesi initializes the array wvhaes which can then be
        !     used repeatedly by subroutine vhaes until nlat or nlon is changed.
        !
        !     input parameters
        !
        !     nlat   the number of colatitudes on the full sphere including the
        !            poles. for example, nlat = 37 for a five degree grid.
        !            nlat determines the grid increment in colatitude as
        !            pi/(nlat-1).  if nlat is odd the equator is located at
        !            grid point i=(nlat+1)/2. if nlat is even the equator is
        !            located half way between points i=nlat/2 and i=nlat/2+1.
        !            nlat must be at least 3. note: on the half sphere, the
        !            number of grid points in the colatitudinal direction is
        !            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than zero. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !     lvhaes the dimension of the array wvhaes as it appears in the
        !            program that calls vhaes. define
        !
        !               l1 = min0(nlat,nlon/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lvhaes must be at least
        !
        !               l1*l2*(nlat+nlat-l1+1)+nlon+15
        !
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls vhaes. lwork must be at least
        !
        !              3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     wvhaes an array which is initialized for use by subroutine vhaes.
        !            once initialized, wvhaes can be used repeatedly by vhaes
        !            as long as nlat or nlon remain unchanged.  wvhaes must not
        !            be altered between calls of vhaes.
        !
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lvhaes
        !            = 4  error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: error_flag
        integer (IP)    :: lvhaes, ldwork
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Compute various workspace dimensions
        !--------------------------------------------------------------------------------

        ldwork = this%Get_ldwork( nlat )
        lvhaes = this%Get_lvhaes( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_VECTOR_ANALYSIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wvhaes ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhaes, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHAES failed in INITIALIZE_VECTOR_ANALYSIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        allocate ( &
            this%dwork(  1:ldwork ), &
            this%wvhaes( 1:lvhaes ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_VECTOR_ANALYSIS_REG'
            write( stderr, '(A)')  trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            wvhaes => this%wvhaes, &
            dwork  => this%dwork, &
            ierror => error_flag &
            )

            call Vhaesi( nlat, nlon, wvhaes, lvhaes, dwork, ldwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_VECTOR_ANALYSIS_REG'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for WVHAES'
                write( stderr, '(A)') 'size( this%WVHAES )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for WORK'
                write( stderr, '(A)') 'size( this%WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_vector_analysis_reg
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_synthesis_reg( this, nlat, nlon )
        !
        !< Purpose:
        !
        ! Sets the various workspace arrays required for
        ! (real) vector harmonic analysis on a regular (equally-spaced) grid
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhsesi.html
        !
        !--------------------------------------------------------------------------------
        !
        !     Documentation: SPHEREPACK 3.2
        !
        !     subroutine vhsesi(nlat,nlon,wvhses,lvhses,work,lwork,ierror)
        !
        !     subroutine vhsesi initializes the array wvhses which can then be
        !     used repeatedly by subroutine vhses until nlat or nlon is changed.
        !
        !     input parameters
        !
        !     nlat   the number of colatitudes on the full sphere including the
        !            poles. for example, nlat = 37 for a five degree grid.
        !            nlat determines the grid increment in colatitude as
        !            pi/(nlat-1).  if nlat is odd the equator is located at
        !            grid point i=(nlat+1)/2. if nlat is even the equator is
        !            located half way between points i=nlat/2 and i=nlat/2+1.
        !            nlat must be at least 3. note: on the half sphere, the
        !            number of grid points in the colatitudinal direction is
        !            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
        !
        !     nlon   the number of distinct londitude points.  nlon determines
        !            the grid increment in longitude as 2*pi/nlon. for example
        !            nlon = 72 for a five degree grid. nlon must be greater
        !            than zero. the axisymmetric case corresponds to nlon=1.
        !            the efficiency of the computation is improved when nlon
        !            is a product of small prime numbers.
        !
        !     lvhses the dimension of the array wvhses as it appears in the
        !            program that calls vhses. define
        !
        !               l1 = min0(nlat,nlon/2) if nlon is even or
        !               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
        !
        !            and
        !
        !               l2 = nlat/2        if nlat is even or
        !               l2 = (nlat+1)/2    if nlat is odd
        !
        !            then lvhses must be at least
        !
        !                  l1*l2*(nlat+nlat-l1+1)+nlon+15
        !
        !
        !     work   a work array that does not have to be saved.
        !
        !     lwork  the dimension of the array work as it appears in the
        !            program that calls vhses. lwork must be at least
        !
        !              3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat
        !
        !
        !     **************************************************************
        !
        !     output parameters
        !
        !     wvhses an array which is initialized for use by subroutine vhses.
        !            once initialized, wvhses can be used repeatedly by vhses
        !            as long as nlat or nlon remain unchanged.  wvhses must not
        !            be altered between calls of vhses.
        !
        !
        !     ierror = 0  no errors
        !            = 1  error in the specification of nlat
        !            = 2  error in the specification of nlon
        !            = 3  error in the specification of lvhses
        !            = 4  error in the specification of lwork
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: error_flag
        integer (IP)    :: ldwork, lvhses
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Compute various workspace dimensions
        !--------------------------------------------------------------------------------

        ldwork = this%Get_ldwork( nlat )
        lvhses = this%Get_lvhses( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Deallocatea arrays ( if necessary )
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%dwork ) ) then

            ! Deallocate array
            deallocate ( &
                this%dwork, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating DWORK failed in INITIALIZE_VECTOR_SYNTHESIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%wvhses ) ) then

            ! Deallocate array
            deallocate ( &
                this%wvhses, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then

                write( stderr, '(A)' ) 'TYPE (workspace_t)'
                write( stderr, '(A)' ) 'Deallocating WVHSES failed in INITIALIZE_VECTOR_SYNTHESIS_REG'
                write( stderr, '(A)' ) trim( error_message )

            end if
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        allocate ( &
            this%dwork(  1:ldwork ), &
            this%wvhses( 1:lvhses ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_VECTOR_SYNTHESIS_REG'
            write( stderr, '(A)')  trim( error_message )

        end if

        !--------------------------------------------------------------------------------
        ! Invoke SPHEREPACK 3.2
        !--------------------------------------------------------------------------------

        associate( &
            wvhses => this%wvhses, &
            dwork  => this%dwork, &
            ierror => error_flag &
            )

            call Vhsesi( nlat, nlon, wvhses, lvhses, dwork, ldwork, ierror )

        end associate

        !--------------------------------------------------------------------------------
        ! Address the error flag
        !--------------------------------------------------------------------------------

        if ( error_flag /= 0 ) then

            write( stderr, '(A)') 'TYPE (workspace_t) in INITIALIZE_VECTOR_SYNTHESIS_REG'

            if ( error_flag == 1 ) then

                write( stderr, '(A)') 'Error in the specification of NLAT'

            else if ( error_flag == 2 ) then

                write( stderr, '(A)') 'Error in the specification of NLON'

            else if ( error_flag == 3 ) then

                write( stderr, '(A)') 'Error in the specification of extent for WVHSES'
                write( stderr, '(A)') 'size( this%WVHSES )'

            else if (error_flag == 4) then

                write( stderr, '(A)') 'Error in the specification of extent for WORK'
                write( stderr, '(A)') 'size( this%WORK )'

            else

                write( stderr, '(A)') 'Undetermined error flag'

            end if
        end if

    end subroutine Initialize_vector_synthesis_reg
    !
    !*****************************************************************************************
    !
    subroutine Initialize_vector_transform_reg( this, nlat, nlon )
        !
        !< Purpose:
        !
        ! Sets the various workspace arrays and pointers
        ! required for the (real) vector harmonic transform
        ! on a regular (equally-spaced) grid
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (workspace_t), intent (in out) :: this
        integer (IP),        intent (in)     :: nlat
        integer (IP),        intent (in)     :: nlon
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check initialization flag
        !--------------------------------------------------------------------------------

        call this%Assert_initialized()

        !--------------------------------------------------------------------------------
        ! Set up vector analysis - REG
        !--------------------------------------------------------------------------------

        call this%Initialize_vector_analysis_reg( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Set up vector synthesis - REG
        !--------------------------------------------------------------------------------

        call this%Initialize_vector_synthesis_reg( nlat, nlon )

        !--------------------------------------------------------------------------------
        ! Allocate pointers for the vector transform coefficients - REG
        !--------------------------------------------------------------------------------

        allocate ( &
            this%real_polar_harmonic_coefficients(          1:nlat, 1:nlat ), &
            this%imaginary_polar_harmonic_coefficients(     1:nlat, 1:nlat ), &
            this%real_azimuthal_harmonic_coefficients(      1:nlat, 1:nlat ), &
            this%imaginary_azimuthal_harmonic_coefficients( 1:nlat, 1:nlat ), &
            stat   = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then

            write( stderr, '(A)' ) 'TYPE (workspace_t)'
            write( stderr, '(A)' ) 'Allocation failed in INITIALIZE_VECTOR_TRANSFORM_REG'
            write( stderr, '(A)' ) trim( error_message )

        end if

    end subroutine Initialize_vector_transform_reg
    !
    !*****************************************************************************************
    !
    pure function Get_lwork( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------

        return_value = (4 * nlon + 2) * nlat

    end function Get_lwork
    !
    !*****************************************************************************************
    !
    pure function Get_ldwork( nlat ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        !--------------------------------------------------------------------------------

        return_value = (3 * nlat * (nlat + 3) + 2)/2

    end function Get_ldwork
    !
    !*****************************************************************************************
    !
    pure function Get_lshags( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            nlat * (3 * (l1 + l2) - 2) &
            +(l1 - 1) * (l2 * (2 * nlat - l1) - 3 * l1)/2 &
            + nlon + 15

    end function Get_lshags
    !
    !*****************************************************************************************
    !
    pure function Get_lshaes( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            (l1 * l2 * (nlat + nlat - l1 + 1))/2 + nlon + 15

    end function Get_lshaes
    !
    !*****************************************************************************************
    !
    pure function Get_lshsgs( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            nlat * (3 * (l1 + l2) - 2) &
            +(l1 - 1) * (l2 * (2 * nlat - l1) - 3 * l1)/2 &
            + nlon + 15

    end function Get_lshsgs
    !
    !*****************************************************************************************
    !
    pure function Get_lshses( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            (l1 * l2 * (nlat + nlat - l1 + 1))/2 + nlon + 15

    end function Get_lshses
    !
    !*****************************************************************************************
    !
    pure function Get_lvhags( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------

        return_value = &
            (nlat + 1) * (nlat + 1) * nlat / 2 &
            + nlon + 15

    end function Get_lvhags
    !
    !*****************************************************************************************
    !
    pure function Get_lvhaes( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        integer (IP) ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            l1 * l2 * (nlat + nlat - l1 + 1) + nlon + 15

    end function Get_lvhaes
    !
    !*****************************************************************************************
    !
    pure function Get_lvhsgs( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)              :: return_value
        integer (IP), intent (in) :: nlat
        integer (IP), intent (in) :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, nlon/2 )
        else
            l1 = min( nlat,(nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value =  &
            l1 * l2 * (nlat + nlat - l1 + 1) &
            + nlon + 15 + 2 * nlat

    end function Get_lvhsgs
    !
    !*****************************************************************************************
    !
    pure function Get_lvhses( nlat, nlon ) result ( return_value )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP)               :: return_value
        integer (IP), intent (in)  :: nlat
        integer (IP), intent (in)  :: nlon
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP) ::  l1, l2
        !--------------------------------------------------------------------------------

        ! Compute parity
        if ( mod( nlon, 2 ) == 0 ) then
            l1 = min( nlat, nlon/2 )
        else
            l1 = min( nlat,(nlon + 1)/2 )
        end if

        if ( mod( nlat, 2 ) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value =  &
            l1 * l2 * (nlat + nlat - l1 + 1) + nlon + 15

    end function Get_lvhses
    !
    !*****************************************************************************************
    !
    subroutine Finalize_workspace( this )
        !
        ! Purpose:
        !< Finalize object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (workspace_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%Destroy()

    end subroutine Finalize_workspace
    !
    !*****************************************************************************************
    !
end module type_workspace_mod
