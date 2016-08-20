module type_GaussianWorkspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Workspace, only: &
        Workspace

    use scalar_analysis_routines, only: &
        ShagsAux

    use module_shsgs, only: &
        ShsgsAux

    use module_vhags, only: &
        VhagsAux

    use module_vhsgs, only: &
        VhsgsAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GaussianWorkspace



    ! Declare derived data type
    type, extends (Workspace), public :: GaussianWorkspace
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
        !----------------------------------------------------------------------
        procedure,         public  :: create => create_gaussian_workspace
        procedure,         public  :: destroy => destroy_gaussian_workspace
        procedure,         private :: initialize_gaussian_scalar_analysis
        procedure,         private :: initialize_gaussian_scalar_synthesis
        procedure,         private :: initialize_gaussian_scalar_transform
        procedure,         private :: initialize_gaussian_vector_analysis
        procedure,         private :: initialize_gaussian_vector_synthesis
        procedure,         private :: initialize_gaussian_vector_transform
        generic,           public  :: assignment (=) => copy_gaussian_workspace
        procedure,         private :: copy_gaussian_workspace
        final                      :: finalize_gaussian_workspace
        !----------------------------------------------------------------------
    end type GaussianWorkspace



    ! Declare constructor
    interface GaussianWorkspace
        module procedure gaussian_workspace_constructor
    end interface



contains



    function gaussian_workspace_constructor(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip),         intent (in) :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip),         intent (in) :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        type (GaussianWorkspace)          :: return_value
        !----------------------------------------------------------------------

        call return_value%create(nlat, nlon)

    end function gaussian_workspace_constructor



    subroutine copy_gaussian_workspace(this, object_to_be_copied)
        !--------------------------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------------------------
        class (GaussianWorkspace), intent (out) :: this
        class (GaussianWorkspace), intent (in)  :: object_to_be_copied
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (object_to_be_copied%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (GaussianWorkspace): '&
                //'in assignment (=) '
        end if

        !
        !==> Make copies
        !
        this%initialized = object_to_be_copied%initialized
        this%legendre_workspace = object_to_be_copied%legendre_workspace
        this%forward_scalar = object_to_be_copied%forward_scalar
        this%forward_vector = object_to_be_copied%forward_vector
        this%backward_scalar = object_to_be_copied%backward_scalar
        this%backward_vector = object_to_be_copied%backward_vector
        this%real_harmonic_coefficients = object_to_be_copied%real_harmonic_coefficients
        this%imaginary_harmonic_coefficients = object_to_be_copied%imaginary_harmonic_coefficients
        this%real_polar_harmonic_coefficients = object_to_be_copied%real_polar_harmonic_coefficients
        this%imaginary_polar_harmonic_coefficients = object_to_be_copied%imaginary_polar_harmonic_coefficients
        this%real_azimuthal_harmonic_coefficients = object_to_be_copied%real_azimuthal_harmonic_coefficients
        this%imaginary_azimuthal_harmonic_coefficients = object_to_be_copied%imaginary_azimuthal_harmonic_coefficients


    end subroutine copy_gaussian_workspace



    subroutine create_gaussian_workspace(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        call this%initialize_gaussian_scalar_transform(nlat, nlon)
        call this%initialize_gaussian_vector_transform(nlat, nlon)
        call get_legendre_workspace(nlat, nlon, this%legendre_workspace)

        ! Set flag
        this%initialized = .true.

    end subroutine create_gaussian_workspace



    subroutine destroy_gaussian_workspace(this)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if (.not.this%initialized) return

        ! Release memory from parent type
        call this%destroy_workspace()

        ! Set flag
        this%initialized = .true.

    end subroutine destroy_gaussian_workspace



    subroutine initialize_gaussian_scalar_analysis(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lshags
        real (wp), allocatable :: work(:), dwork(:)
        type (ShagsAux)        :: aux
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        lwork = get_lwork(nlat)
        ldwork = get_ldwork(nlat)
        lshags = aux%get_lshags(nlat, nlon)

        !
        !==> Allocate memory
        !
        if (allocated(this%forward_scalar)) deallocate(this%forward_scalar )

        !
        !==> Allocate memory
        !
        allocate( this%forward_scalar(lshags) )
        allocate( work(lwork) )
        allocate( dwork(ldwork) )

        associate( &
            wshags => this%forward_scalar, &
            ierror => error_flag &
            )
            !
            !==> Initialize workspace array for analysis
            !
            call aux%shagsi(nlat, nlon, wshags, lshags, work, lwork, dwork, ldwork, ierror)

        end associate


        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for legendre_workspace'
            case(5)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for dwork'
            case(6)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in call to compute_gaussian_latitudes_and_weights to compute gaussian points '&
                    //'due to failure in eigenvalue routine'
            case default
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'Undetermined error flag'
        end select


        !
        !==> Release memory
        !
        deallocate( work )
        deallocate( dwork )

    end subroutine initialize_gaussian_scalar_analysis



    subroutine initialize_gaussian_scalar_synthesis(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lshsgs
        real (wp), allocatable :: work(:), dwork(:)
        type (ShsgsAux)        :: aux
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        lwork = get_lwork(nlat)
        ldwork = get_ldwork(nlat)
        lshsgs = aux%get_lshsgs(nlat, nlon)

        !
        !==> Allocate memory
        !
        if (allocated(this%backward_scalar)) deallocate(this%backward_scalar )

        allocate( work(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%backward_scalar(lshsgs) )


        associate( &
            wshsgs => this%backward_scalar, &
            ierror => error_flag &
            )
            !
            !==> Initialize workspace array for synthesis
            !
            call aux%shsgsi(nlat, nlon, wshsgs, lshsgs, work, lwork, dwork, ldwork, ierror)

        end associate

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for backward_scalar'
            case(4)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for legendre_workspace'
            case(5)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for dwork'
            case(6)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in call to compute_gaussian_latitudes_and_weights due to failure in eigenvalue routine'
            case default
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'Undetermined error flag'
        end select

        !
        !==> Release memory
        !
        deallocate( work )
        deallocate( dwork )

    end subroutine initialize_gaussian_scalar_synthesis



    subroutine initialize_gaussian_scalar_transform(this, nlat, nlon)
        !
        ! Purpose:
        !
        !  Set the various workspace arrays and pointers to perform the
        !  (real) scalar harmonic transform on a gaussian grid
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Set up scalar analysis
        call this%initialize_gaussian_scalar_analysis(nlat, nlon)

        ! Set up scalar synthesis
        call this%initialize_gaussian_scalar_synthesis(nlat, nlon)

        !
        !==> Allocate memory
        !
        allocate( this%real_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_harmonic_coefficients(nlat, nlat) )

    end subroutine initialize_gaussian_scalar_transform



    subroutine initialize_gaussian_vector_analysis(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: ldwork, lvhags
        real (wp), allocatable :: dwork(:)
        type (VhagsAux)        :: aux
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        ldwork = get_ldwork(nlat)
        lvhags = aux%get_lvhags(nlat, nlon)

        !
        !==> Allocate memory
        !
        if (allocated(this%forward_vector)) deallocate( this%forward_vector )

        allocate( dwork(ldwork) )
        allocate(this%forward_vector(lvhags) )


        associate( &
            wvhags => this%forward_vector, &
            ierror => error_flag &
            )

            !
            !==> Initialize workspace arrays for vector analysis
            !
            call aux%vhagsi(nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror)

        end associate


        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of extent for forward_vector'
            case(4)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'error in the specification of extent for dwork'
            case default
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_analysis'&
                    //'Undetermined error flag'
        end select

        !
        !==> Release memory
        !
        deallocate( dwork )

    end subroutine initialize_gaussian_vector_analysis



    subroutine initialize_gaussian_vector_synthesis(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: ldwork, lvhsgs
        real (wp), allocatable :: dwork(:)
        type (VhsgsAux)        :: aux
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        ldwork = get_ldwork(nlat)
        lvhsgs = aux%get_lvhsgs(nlat, nlon)

        !
        !==> Allocate memory
        !
        if (allocated(this%backward_vector)) deallocate( this%backward_vector )

        allocate( dwork(ldwork) )
        allocate( this%backward_vector(lvhsgs) )


        associate( &
            wvhsgs => this%backward_vector, &
            ierror => error_flag &
            )
            !
            !==> Initialize workspace arrays for vector synthesis
            !
            call aux%vhsgsi(nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror)

        end associate

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of extent for backward_vector'
            case(4)
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of extent for dwork'
            case default
                error stop 'Object of class (GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'Undetermined error flag'
        end select

        !
        !==> Release memory
        !
        deallocate( dwork )

    end subroutine initialize_gaussian_vector_synthesis



    subroutine initialize_gaussian_vector_transform(this, nlat, nlon)
        !
        ! Purpose:
        !
        !  Sets the various workspace arrays and pointers
        !  required for the (real) vector harmonic transform
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Set up vector analysis
        call this%initialize_gaussian_vector_analysis(nlat, nlon)

        ! Set up vector analysis
        call this%initialize_gaussian_vector_synthesis(nlat, nlon)

        !
        !==> Allocate memory
        !
        allocate( this%real_polar_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_polar_harmonic_coefficients(nlat, nlat) )
        allocate( this%real_azimuthal_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_azimuthal_harmonic_coefficients(nlat, nlat) )

    end subroutine initialize_gaussian_vector_transform



    pure function get_lwork(nlat) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat
        integer (ip)              :: return_value
        !----------------------------------------------------------------------
        type (ShagsAux) :: shags_aux
        type (ShsgsAux) :: shsgs_aux
        integer (ip)    :: lwork(2)
        !----------------------------------------------------------------------

        lwork(1) = shags_aux%get_lwork(nlat)
        lwork(2) = shsgs_aux%get_lwork(nlat)

        return_value = maxval(lwork)

    end function get_lwork



    pure function get_ldwork(nlat) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        type (ShagsAux) :: shags_aux
        type (ShsgsAux) :: shsgs_aux
        type (VhagsAux) :: vhags_aux
        type (VhsgsAux) :: vhsgs_aux
        integer (ip)    :: ldwork(4)
        !----------------------------------------------------------------------

        ldwork(1) = shags_aux%get_ldwork(nlat)
        ldwork(2) = shsgs_aux%get_ldwork(nlat)
        ldwork(3) = vhags_aux%get_ldwork(nlat)
        ldwork(4) = vhsgs_aux%get_ldwork(nlat)

        return_value = maxval(ldwork)

    end function get_ldwork


    pure subroutine get_legendre_workspace(nlat, nlon, workspace, nt, ityp, isym)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip),           intent (in)  :: nlat
        integer (ip),           intent (in)  :: nlon
        real (wp), allocatable, intent (out) :: workspace(:)
        integer (ip), optional, intent (in)  :: nt
        integer (ip), optional, intent (in)  :: ityp
        integer (ip), optional, intent (in)  :: isym
        !----------------------------------------------------------------------
        type (ShagsAux) :: shags_aux
        type (ShsgsAux) :: shsgs_aux
        type (VhagsAux) :: vhags_aux
        type (VhsgsAux) :: vhsgs_aux
        integer (ip)    :: work_size(4)
        integer (ip)    :: lwork, nt_op, ityp_op, isym_op
        !----------------------------------------------------------------------

        !
        !==> Address optional arguments
        !
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        if (present(isym)) then
            isym_op = isym
        else
            isym_op = 0
        end if

        if (present(ityp)) then
            ityp_op = ityp
        else
            ityp_op = 0
        end if

        work_size(1) = shags_aux%get_legendre_workspace_size(nlat, nlon, nt_op, isym_op)
        work_size(2) = shsgs_aux%get_legendre_workspace_size(nlat, nlon, nt_op, isym_op)
        work_size(3) = vhags_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(4) = vhsgs_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)

        lwork = maxval(work_size)
        !
        !==> Allocate memory
        !
        allocate( workspace(lwork) )

    end subroutine get_legendre_workspace



    subroutine finalize_gaussian_workspace(this)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        type (GaussianWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_gaussian_workspace



end module type_GaussianWorkspace
