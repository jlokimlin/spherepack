module type_RegularWorkspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Workspace, only: &
        Workspace

    use module_shaes, only: &
        ShaesAux

    use module_shses, only: &
        ShsesAux

    use module_vhaes, only: &
        VhaesAux

    use module_vhses, only: &
        VhsesAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: RegularWorkspace


    ! Declare derived data type
    type, extends (Workspace), public :: RegularWorkspace
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public  :: create => create_regular_workspace
        procedure, public  :: destroy => destroy_regular_workspace
        procedure, private :: initialize_regular_scalar_analysis
        procedure, private :: initialize_regular_scalar_synthesis
        procedure, private :: initialize_regular_scalar_transform
        procedure, private :: initialize_regular_vector_analysis
        procedure, private :: initialize_regular_vector_synthesis
        procedure, private :: initialize_regular_vector_transform
        generic,   public  :: assignment (=) => copy_regular_workspace
        procedure, private :: copy_regular_workspace
        final              :: finalize_regular_workspace
        !----------------------------------------------------------------------
    end type RegularWorkspace



    ! Declare constructor
    interface RegularWorkspace
        module procedure regular_workspace_constructor
    end interface



contains



    function regular_workspace_constructor(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip), intent (in) :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        type (RegularWorkspace)   :: return_value
        !----------------------------------------------------------------------

        call return_value%create(nlat, nlon)

    end function regular_workspace_constructor



    subroutine copy_regular_workspace(this, object_to_be_copied)
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (RegularWorkspace), intent (out) :: this
        class (RegularWorkspace), intent (in)  :: object_to_be_copied
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if (object_to_be_copied%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (RegularWorkspace): '&
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


    end subroutine copy_regular_workspace



    subroutine create_regular_workspace(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        integer (ip),             intent (in)     :: nlat
        integer (ip),             intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set up transforms for regular grids
        call this%initialize_regular_scalar_transform(nlat, nlon)
        call this%initialize_regular_vector_transform(nlat, nlon)
        call get_legendre_workspace(nlat, nlon, this%legendre_workspace)

        ! Set flag
        this%initialized = .true.

    end subroutine create_regular_workspace



    subroutine destroy_regular_workspace(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if (.not.this%initialized) return

        ! Release memory from parent type
        call this%destroy_workspace()

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_regular_workspace



    subroutine initialize_regular_scalar_analysis(this, nlat, nlon)
        !
        ! Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a regular (equally-spaced) grid.
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shaesi.html
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        integer (ip),             intent (in)     :: nlat
        integer (ip),             intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lshaes
        real (wp), allocatable :: dwork(:), work(:)
        type (ShaesAux)        :: aux
        !----------------------------------------------------------------------

        !
        !==> Compute dimensions of various workspace arrays
        !
        lwork = get_lwork(nlat, nlon)
        ldwork = get_ldwork(nlat)
        lshaes = aux%get_lshaes(nlat, nlon)

        !
        !==>  Allocate memory
        !
        if (allocated(this%forward_scalar)) deallocate( this%forward_scalar )
        allocate( work(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%forward_scalar(lshaes) )


        associate( &
            wshaes => this%forward_scalar, &
            ierror => error_flag &
            )
            !
            !==> Initialize workspace for scalar synthesis
            !
            call aux%shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, ldwork, ierror)

        end associate

        !
        !==>  Address error flag
        !
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of extent for legendre_workspace'
            case default
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'Undetermined error flag'
        end select

        !
        !==> Release memory
        !
        deallocate( work )
        deallocate( dwork )

    end subroutine initialize_regular_scalar_analysis



    subroutine initialize_regular_scalar_synthesis(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        integer (ip),             intent (in)     :: nlat
        integer (ip),             intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lshses
        real (wp), allocatable :: dwork(:), work(:)
        type (ShsesAux)        :: aux
        !----------------------------------------------------------------------

        ! Set up various workspace dimensions
        lwork = get_lwork(nlat, nlon)
        ldwork = get_ldwork(nlat)
        lshses = aux%get_lshses(nlat, nlon)

        !
        !==> Allocate memory
        !
        if (allocated(this%backward_scalar)) deallocate( this%backward_scalar )
        allocate( this%backward_scalar(lshses) )
        allocate( work(lwork) )
        allocate( dwork(ldwork) )


        associate( &
            wshses => this%backward_scalar, &
            ierror => error_flag &
            )
            !
            !==> Initialize workspace for scalar synthesis
            !
            call aux%shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, ldwork, ierror)

        end associate


        !  Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of extent for legendre_workspace'
            case default
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'Undetermined error flag'
        end select

        !
        !==> Release memory
        !
        deallocate( work )
        deallocate( dwork )

    end subroutine initialize_regular_scalar_synthesis



    subroutine initialize_regular_scalar_transform(this, nlat, nlon)
        !
        ! Purpose:
        !
        !  Set the various workspace arrays and pointers to perform the
        !  (real) scalar harmonic transform on a regular (equally-spaced) grid.
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        integer (ip),             intent (in)     :: nlat
        integer (ip),             intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Set up scalar analysis
        call this%initialize_regular_scalar_analysis(nlat, nlon)

        ! Set up scalar synthesis
        call this%initialize_regular_scalar_synthesis(nlat, nlon)

        ! Allocate memory for the (real) scalar harmonic transform
        allocate( this%real_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_harmonic_coefficients(nlat, nlat) )

    end subroutine initialize_regular_scalar_transform



    subroutine initialize_regular_vector_analysis(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        integer (ip),             intent (in)     :: nlat
        integer (ip),             intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lvhaes
        real (wp), allocatable :: work(:), dwork(:)
        type (VhaesAux)        :: aux
        !----------------------------------------------------------------------

        ! Compute various workspace dimensions
        lwork = get_lwork(nlat, nlon)
        ldwork = get_ldwork(nlat)
        lvhaes = aux%get_lvhaes(nlat, nlon)

        !
        !==> Allocate memory
        !
        if (allocated(this%forward_vector)) deallocate( this%forward_vector )
        allocate( work(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%forward_vector(lvhaes) )


        associate( &
            wvhaes => this%forward_vector, &
            ierror => error_flag &
            )
            !
            !==> Initialize workspace for analysis
            !
            call aux%vhaesi(nlat, nlon, wvhaes, lvhaes, work, lwork, dwork, ldwork, ierror)

        end associate

        ! Address the error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of extent for forward_vector'
            case(4)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of extent for unsaved work'
            case(5)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //'error in the specification of extent for unsaved dwork'
            case default
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_analysis'&
                    //' Undetermined error flag'
        end select

        !
        !==> Release memory
        !
        deallocate( work )
        deallocate( dwork )

    end subroutine initialize_regular_vector_analysis



    subroutine initialize_regular_vector_synthesis(this, nlat, nlon)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        integer (ip),             intent (in)     :: nlat
        integer (ip),             intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lvhses
        real (wp), allocatable :: work(:), dwork(:)
        type (VhsesAux)        :: aux
        !----------------------------------------------------------------------

        ! Compute various workspace dimensions
        lwork = get_lwork(nlat, nlon)
        ldwork = get_ldwork(nlat)
        lvhses = aux%get_lvhses(nlat, nlon)

        !
        !==> Allocate memory
        !
        if (allocated(this%backward_vector)) deallocate( this%backward_vector )
        allocate( work(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%backward_vector(lvhses) )

        associate( &
            wvhses => this%backward_vector, &
            ierror => error_flag &
            )

            !
            !==> Initialize workspace for vector synthesis
            !
            call aux%vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, ldwork, ierror)

        end associate

        ! Address the error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of extent for backward_vector'
            case(4)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis'&
                    //'error in the specification of extent for unsaved work'
            case(5)
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis '&
                    //'error in the specification of extent for unsaved dwork'
            case default
                error stop 'Object of class (RegularWorkspace): '&
                    //'in initialize_regular_vector_synthesis'&
                    //' Undetermined error flag'
        end select

        !
        !==> Release memory
        !
        deallocate( work )
        deallocate( dwork )

    end subroutine initialize_regular_vector_synthesis



    subroutine initialize_regular_vector_transform(this, nlat, nlon)
        !
        ! Purpose:
        !
        ! Sets the various workspace arrays and pointers
        ! required for the (real) vector harmonic transform
        ! on a regular (equally-spaced) grid
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        integer (ip),             intent (in)     :: nlat
        integer (ip),             intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Set up vector analysis
        call this%initialize_regular_vector_analysis(nlat, nlon)

        ! Set up vector synthesis
        call this%initialize_regular_vector_synthesis(nlat, nlon)

        ! Allocate memory for the vector transform coefficients
        allocate( this%real_polar_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_polar_harmonic_coefficients(nlat, nlat) )
        allocate( this%real_azimuthal_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_azimuthal_harmonic_coefficients(nlat, nlat) )

    end subroutine initialize_regular_vector_transform



    pure function get_lwork(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        type (ShaesAux) :: shaes_aux
        type (ShsesAux) :: shses_aux
        type (VhaesAux) :: vhaes_aux
        type (VhsesAux) :: vhses_aux
        integer (ip)    :: lwork(4)
        !----------------------------------------------------------------------

        lwork(1) = shaes_aux%get_lwork(nlat, nlon)
        lwork(2) = shses_aux%get_lwork(nlat, nlon)
        lwork(3) = vhaes_aux%get_lwork(nlat, nlon)
        lwork(4) = vhses_aux%get_lwork(nlat, nlon)

        return_value = maxval(lwork)

    end function get_lwork


    pure function get_ldwork(nlat) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        type (ShaesAux) :: shaes_aux
        type (ShsesAux) :: shses_aux
        type (VhaesAux) :: vhaes_aux
        type (VhsesAux) :: vhses_aux
        integer (ip)    :: ldwork(4)
        !----------------------------------------------------------------------

        ldwork(1) = shaes_aux%get_ldwork(nlat)
        ldwork(2) = shses_aux%get_ldwork(nlat)
        ldwork(3) = vhaes_aux%get_ldwork(nlat)
        ldwork(4) = vhses_aux%get_ldwork(nlat)

        return_value = maxval(ldwork)

    end function get_ldwork


    pure subroutine get_legendre_workspace(nlat, nlon, workspace, nt, ityp)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip),           intent (in)  :: nlat
        integer (ip),           intent (in)  :: nlon
        real (wp), allocatable, intent (out) :: workspace(:)
        integer (ip), optional, intent (in)  :: nt
        integer (ip), optional, intent (in)  :: ityp
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        type (ShaesAux) :: shaes_aux
        type (ShsesAux) :: shses_aux
        type (VhaesAux) :: vhaes_aux
        type (VhsesAux) :: vhses_aux
        integer (ip)    :: work_size(4)
        integer (ip)    :: lwork, nt_op, ityp_op
        !----------------------------------------------------------------------

        !
        !==> Address optional arguments
        !
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        if (present(ityp)) then
            ityp_op = ityp
        else
            ityp_op = 0
        end if

        work_size(1) = shaes_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(2) = shses_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(3) = vhaes_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(4) = vhses_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)

        lwork = maxval(work_size)
        !
        !==> Allocate memory
        !
        allocate( workspace(lwork) )

    end subroutine get_legendre_workspace



    subroutine finalize_regular_workspace(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (RegularWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy_workspace()

    end subroutine finalize_regular_workspace



end module type_RegularWorkspace
