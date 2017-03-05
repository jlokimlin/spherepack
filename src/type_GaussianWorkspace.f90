module type_GaussianWorkspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Workspace, only: &
        Workspace

    use scalar_analysis_routines, only: &
        ShagsAux

    use scalar_synthesis_routines, only: &
        ShsgsAux

    use vector_analysis_routines, only: &
        VhagsAux

    use vector_synthesis_routines, only: &
        VhsgsAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GaussianWorkspace

    type, public, extends(Workspace) :: GaussianWorkspace
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_gaussian_workspace
        procedure, public  :: destroy => destroy_gaussian_workspace
        procedure, private :: initialize_gaussian_scalar_analysis
        procedure, private :: initialize_gaussian_scalar_synthesis
        procedure, private :: initialize_gaussian_vector_analysis
        procedure, private :: initialize_gaussian_vector_synthesis
        procedure, private :: copy_gaussian_workspace
        ! Generic type-bound procedures
        generic, public :: assignment (=) => copy_gaussian_workspace
    end type GaussianWorkspace

    ! Declare user-defined constructor
    interface GaussianWorkspace
        module procedure gaussian_workspace_constructor
    end interface

contains

    function gaussian_workspace_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(GaussianWorkspace)           :: return_value

        ! Local variables
        integer(ip) :: nt_op

        ! Address optional argument
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        call return_value%create(nlat, nlon, nt_op)

    end function gaussian_workspace_constructor

    subroutine copy_gaussian_workspace(self, other)

        ! Dummy arguments
        class(GaussianWorkspace), intent(out) :: self
        class(GaussianWorkspace), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(GaussianWorkspace): '&
                //'in assignment (=) '
        end if

        !  Make copies
        call self%copy_workspace(other)

    end subroutine copy_gaussian_workspace

    subroutine create_gaussian_workspace(self, nlat, nlon, nt)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon
        integer(ip),              intent(in)    :: nt

        ! Ensure that object is usable
        call self%destroy()

        ! Initialize harmonic coefficients
        call self%initialize_harmonic_coefficients(nlat, nlon, nt)

        ! Set up scalar analysis
        call self%initialize_gaussian_scalar_analysis(nlat, nlon)

        ! Set up scalar synthesis
        call self%initialize_gaussian_scalar_synthesis(nlat, nlon)

        ! Set up vector analysis
        call self%initialize_gaussian_vector_analysis(nlat, nlon)

        ! Set up vector analysis
        call self%initialize_gaussian_vector_synthesis(nlat, nlon)

        call get_legendre_workspace(nlat, nlon, self%legendre_workspace)

        ! Set flag
        self%initialized = .true.

    end subroutine create_gaussian_workspace

    subroutine destroy_gaussian_workspace(self)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        ! Release memory from parent type
        call self%destroy_workspace()

        ! Set flag
        self%initialized = .true.

    end subroutine destroy_gaussian_workspace

    subroutine initialize_gaussian_scalar_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon

        ! Local variables
        integer(ip)    :: error_flag, lshags
        type(ShagsAux) :: aux

        ! Compute dimensions of various workspace arrays
        lshags = aux%get_lshags(nlat, nlon)

        !  Allocate memory
        if (allocated(self%forward_scalar)) deallocate (self%forward_scalar)
        allocate (self%forward_scalar(lshags))

        call aux%shagsi(nlat, nlon, self%forward_scalar, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for legendre_workspace'
            case(5)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in the specification of extent for dwork'
            case(6)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'error in call to compute_gaussian_latitudes_and_weights to compute gaussian points '&
                    //'due to failure in eigenvalue routine'
            case default
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_analysis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_scalar_analysis

    subroutine initialize_gaussian_scalar_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon

        ! Local variables
        integer(ip)    :: error_flag
        integer(ip)    :: lshsgs
        type(ShsgsAux) :: aux

        ! Compute dimensions of various workspace arrays
        lshsgs = aux%get_lshsgs(nlat, nlon)

        !  Allocate memory
        if (allocated(self%backward_scalar)) deallocate (self%backward_scalar)
        allocate (self%backward_scalar(lshsgs))

        call aux%shsgsi(nlat, nlon, self%backward_scalar, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for backward_scalar'
            case(4)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for legendre_workspace'
            case(5)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in the specification of extent for dwork'
            case(6)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'error in call to compute_gaussian_latitudes_and_weights due to failure in eigenvalue routine'
            case default
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_scalar_synthesis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_scalar_synthesis

    subroutine initialize_gaussian_vector_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout)  :: self
        integer(ip),              intent(in)     :: nlat
        integer(ip),              intent(in)     :: nlon

        ! Local variables
        integer(ip)     :: error_flag, lvhags
        type(VhagsAux)  :: aux

        ! Compute dimensions of various workspace arrays
        lvhags = aux%get_lvhags(nlat, nlon)

        !  Allocate memory
        if (allocated(self%forward_vector)) deallocate (self%forward_vector)
        allocate(self%forward_vector(lvhags))

            call aux%vhagsi(nlat, nlon, self%forward_vector, error_flag)

            ! Address error flag
            select case (error_flag)
                case(0)
                    return
                case(1)
                    error stop 'Object of class(GaussianWorkspace) '&
                        //'in initialize_gaussian_vector_analysis'&
                        //'error in the specification of NUMBER_OF_LATITUDES'
                case(2)
                    error stop 'Object of class(GaussianWorkspace) '&
                        //'in initialize_gaussian_vector_analysis'&
                        //'error in the specification of NUMBER_OF_LONGITUDES'
                case(3)
                    error stop 'Object of class(GaussianWorkspace) '&
                        //'in initialize_gaussian_vector_analysis'&
                        //'error in the specification of extent for forward_vector'
                case(4)
                    error stop 'Object of class(GaussianWorkspace) '&
                        //'in initialize_gaussian_vector_analysis'&
                        //'error in the specification of extent for dwork'
                case default
                    error stop 'Object of class(GaussianWorkspace) '&
                        //'in initialize_gaussian_vector_analysis'&
                        //'Undetermined error flag'
            end select

    end subroutine initialize_gaussian_vector_analysis

    subroutine initialize_gaussian_vector_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(GaussianWorkspace), intent(inout) :: self
        integer(ip),              intent(in)    :: nlat
        integer(ip),              intent(in)    :: nlon

        ! Local variables
        integer(ip)    :: error_flag, lvhsgs
        type(VhsgsAux) :: aux

        ! Compute required workspace sizes
        lvhsgs = aux%get_lvhsgs(nlat, nlon)

        !  Allocate memory
        if (allocated(self%backward_vector)) deallocate (self%backward_vector)
        allocate (self%backward_vector(lvhsgs))

        call aux%vhsgsi(nlat, nlon, self%backward_vector, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of extent for backward_vector'
            case(4)
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'error in the specification of extent for dwork'
            case default
                error stop 'Object of class(GaussianWorkspace) '&
                    //'in initialize_gaussian_vector_synthesis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_vector_synthesis

    pure function get_lwork(nlat) result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip)             :: return_value

        ! Local variables
        type(ShagsAux) :: shags_aux
        type(ShsgsAux) :: shsgs_aux
        integer(ip)    :: lwork(2)

        lwork(1) = shags_aux%get_lwork(nlat)
        lwork(2) = shsgs_aux%get_lwork(nlat)

        return_value = maxval(lwork)

    end function get_lwork

    pure function get_ldwork(nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip)             :: return_value

        ! Local variables
        type(ShagsAux) :: shags_aux
        type(ShsgsAux) :: shsgs_aux
        type(VhagsAux) :: vhags_aux
        type(VhsgsAux) :: vhsgs_aux
        integer(ip)    :: ldwork(4)

        ldwork(1) = shags_aux%get_ldwork(nlat)
        ldwork(2) = shsgs_aux%get_ldwork(nlat)
        ldwork(3) = vhags_aux%get_ldwork(nlat)
        ldwork(4) = vhsgs_aux%get_ldwork(nlat)

        return_value = maxval(ldwork)

    end function get_ldwork

    pure subroutine get_legendre_workspace(nlat, nlon, workspace, nt, ityp, isym)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: workspace(:)
        integer(ip), optional, intent(in)  :: nt
        integer(ip), optional, intent(in)  :: ityp
        integer(ip), optional, intent(in)  :: isym

        ! Local variables
        type(ShagsAux) :: shags_aux
        type(ShsgsAux) :: shsgs_aux
        type(VhagsAux) :: vhags_aux
        type(VhsgsAux) :: vhsgs_aux
        integer(ip)    :: work_size(4)
        integer(ip)    :: lwork, nt_op, ityp_op, isym_op

        !  Address optional arguments
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

        ! Compute required workspace size
        work_size(1) = shags_aux%get_legendre_workspace_size(nlat, nlon, nt_op, isym_op)
        work_size(2) = shsgs_aux%get_legendre_workspace_size(nlat, nlon, nt_op, isym_op)
        work_size(3) = vhags_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(4) = vhsgs_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        lwork = maxval(work_size)

        !  Allocate memory
        allocate (workspace(lwork))

    end subroutine get_legendre_workspace

end module type_GaussianWorkspace
