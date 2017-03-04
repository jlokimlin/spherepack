module type_RegularWorkspace

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_Workspace, only: &
        Workspace

    use scalar_analysis_routines, only: &
        ShaesAux

    use scalar_synthesis_routines, only: &
        ShsesAux

    use vector_analysis_routines, only: &
        VhaesAux

    use vector_synthesis_routines, only: &
        VhsesAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: RegularWorkspace
    
    type, public, extends(Workspace) :: RegularWorkspace
    contains
        ! Type-bound procedures
        procedure, public  :: create => create_regular_workspace
        procedure, public  :: destroy => destroy_regular_workspace
        procedure, private :: initialize_regular_scalar_analysis
        procedure, private :: initialize_regular_scalar_synthesis
        procedure, private :: initialize_regular_vector_analysis
        procedure, private :: initialize_regular_vector_synthesis
        procedure, private :: copy_regular_workspace
        ! Generic type-bound procedures
        generic, public :: assignment (=) => copy_regular_workspace
    end type RegularWorkspace

    ! Declare user-defined constructor
    interface RegularWorkspace
        module procedure regular_workspace_constructor
    end interface

contains

    function regular_workspace_constructor(nlat, nlon, nt) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat ! number of latitudinal points 0 <= theta <= pi
        integer(ip),           intent(in) :: nlon ! number of longitudinal points 0 <= phi <= 2*pi
        integer(ip), optional, intent(in) :: nt ! Number of syntheses
        type(RegularWorkspace)            :: return_value

        ! Local variables
        integer(ip) :: nt_op

        ! Address optional argument
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        call return_value%create(nlat, nlon, nt_op)

    end function regular_workspace_constructor

    subroutine copy_regular_workspace(self, other)

        ! Dummy arguments
        class(RegularWorkspace), intent(out) :: self
        class(RegularWorkspace), intent(in)  :: other

        ! Check if object is usable
        if (.not.other%initialized) then
            error stop 'Uninitialized object of class(RegularWorkspace): '&
                //'in assignment (=) '
        end if

        !  Make copies
        call self%copy_workspace(other)

    end subroutine copy_regular_workspace

    subroutine create_regular_workspace(self, nlat, nlon, nt)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon
        integer(ip),             intent(in)    :: nt

        ! Ensure that object is usable
        call self%destroy()

        ! Initialize harmonic coefficients
        call self%initialize_harmonic_coefficients(nlat, nlon, nt)

        ! Set up scalar analysis
        call self%initialize_regular_scalar_analysis(nlat, nlon)

        ! Set up scalar synthesis
        call self%initialize_regular_scalar_synthesis(nlat, nlon)

        ! Set up vector analysis
        call self%initialize_regular_vector_analysis(nlat, nlon)

        ! Set up vector synthesis
        call self%initialize_regular_vector_synthesis(nlat, nlon)

        call get_legendre_workspace(nlat, nlon, self%legendre_workspace)

        ! Set flag
        self%initialized = .true.

    end subroutine create_regular_workspace

    subroutine destroy_regular_workspace(self)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self

        ! Check flag
        if (.not.self%initialized) return

        ! Release memory from parent type
        call self%destroy_workspace()

        ! Reset initialization flag
        self%initialized = .false.

    end subroutine destroy_regular_workspace

    subroutine initialize_regular_scalar_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon

        ! Local variables
        integer(ip)    :: error_flag
        integer(ip)    :: lshaes
        type(ShaesAux) :: aux

        ! Compute dimensions of various workspace arrays
        lshaes = aux%get_lshaes(nlat, nlon)

        ! Allocate memory
        if (allocated(self%forward_scalar)) deallocate( self%forward_scalar )
        allocate( self%forward_scalar(lshaes) )

        ! Call procedural routine
        call aux%shaesi(nlat, nlon, self%forward_scalar, error_flag)

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'error in the specification of extent for legendre_workspace'
            case default
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_analysis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_regular_scalar_analysis

    subroutine initialize_regular_scalar_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout)  :: self
        integer(ip),             intent(in)     :: nlat
        integer(ip),             intent(in)     :: nlon

        ! Local variables
        integer(ip)    :: error_flag, lshses
        type(ShsesAux) :: aux

        ! Set up various workspace dimensions
        lshses = aux%get_lshses(nlat, nlon)

        !  Allocate memory
        if (allocated(self%backward_scalar)) deallocate( self%backward_scalar )
        allocate( self%backward_scalar(lshses) )

        call aux%shsesi(nlat, nlon, self%backward_scalar, error_flag)

        !  Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of extent for forward_scalar'
            case(4)
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'error in the specification of extent for legendre_workspace'
            case default
                error stop 'Object of class(RegularWorkspace): '&
                    //'in initialize_regular_scalar_synthesis '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_regular_scalar_synthesis

    subroutine initialize_regular_vector_analysis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon

        ! Local variables
        integer(ip)     :: error_flag
        integer(ip)     :: lwork, ldwork, lvhaes
        type(VhaesAux)  :: aux


        ! Compute various workspace dimensions
        lwork = get_lwork(nlat, nlon)
        ldwork = get_ldwork(nlat)
        lvhaes = aux%get_lvhaes(nlat, nlon)

        !  Allocate memory
        if (allocated(self%forward_vector)) deallocate( self%forward_vector )
        allocate( self%forward_vector(lvhaes) )

        ! Initialize workspace for analysis
        block
            real(wp) :: work(lwork), dwork(ldwork)

            call aux%vhaesi( &
                nlat, nlon, self%forward_vector, lvhaes, work, lwork, dwork, ldwork, error_flag)

            ! Address the error flag
            select case (error_flag)
                case(0)
                    return
                case(1)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_analysis'&
                        //'error in the specification of NUMBER_OF_LATITUDES'
                case(2)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_analysis'&
                        //'error in the specification of NUMBER_OF_LONGITUDES'
                case(3)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_analysis'&
                        //'error in the specification of extent for forward_vector'
                case(4)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_analysis'&
                        //'error in the specification of extent for unsaved work'
                case(5)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_analysis'&
                        //'error in the specification of extent for unsaved dwork'
                case default
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_analysis'&
                        //' Undetermined error flag'
            end select
        end block

    end subroutine initialize_regular_vector_analysis

    subroutine initialize_regular_vector_synthesis(self, nlat, nlon)

        ! Dummy arguments
        class(RegularWorkspace), intent(inout) :: self
        integer(ip),             intent(in)    :: nlat
        integer(ip),             intent(in)    :: nlon

        ! Local variables
        integer(ip)    :: error_flag
        integer(ip)    :: lwork, ldwork, lvhses
        type(VhsesAux) :: aux

        ! Compute various workspace dimensions
        lwork = get_lwork(nlat, nlon)
        ldwork = get_ldwork(nlat)
        lvhses = aux%get_lvhses(nlat, nlon)

        ! Allocate memory
        if (allocated(self%backward_vector)) deallocate( self%backward_vector )
        allocate( self%backward_vector(lvhses) )

        ! Initialize workspace for vector synthesis
        block
            real(wp) :: work(lwork), dwork(ldwork)

            call aux%vhsesi( &
                nlat, nlon, self%backward_vector, lvhses, work, lwork, dwork, ldwork, error_flag)

            ! Address the error flag
            select case (error_flag)
                case(0)
                    return
                case(1)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_synthesis '&
                        //'error in the specification of NUMBER_OF_LATITUDES'
                case(2)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_synthesis '&
                        //'error in the specification of NUMBER_OF_LONGITUDES'
                case(3)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_synthesis '&
                        //'error in the specification of extent for backward_vector'
                case(4)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_synthesis'&
                        //'error in the specification of extent for unsaved work'
                case(5)
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_synthesis '&
                        //'error in the specification of extent for unsaved dwork'
                case default
                    error stop 'Object of class(RegularWorkspace): '&
                        //'in initialize_regular_vector_synthesis'&
                        //' Undetermined error flag'
            end select
        end block

    end subroutine initialize_regular_vector_synthesis

    pure function get_lwork(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        type(ShaesAux) :: shaes_aux
        type(ShsesAux) :: shses_aux
        type(VhaesAux) :: vhaes_aux
        type(VhsesAux) :: vhses_aux
        integer(ip)    :: lwork(4)

        lwork(1) = shaes_aux%get_lwork(nlat, nlon)
        lwork(2) = shses_aux%get_lwork(nlat, nlon)
        lwork(3) = vhaes_aux%get_lwork(nlat, nlon)
        lwork(4) = vhses_aux%get_lwork(nlat, nlon)

        return_value = maxval(lwork)

    end function get_lwork

    pure function get_ldwork(nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip)             :: return_value

        ! Local variables
        type(ShaesAux) :: shaes_aux
        type(ShsesAux) :: shses_aux
        type(VhaesAux) :: vhaes_aux
        type(VhsesAux) :: vhses_aux
        integer(ip)    :: ldwork(4)

        ldwork(1) = shaes_aux%get_ldwork(nlat)
        ldwork(2) = shses_aux%get_ldwork(nlat)
        ldwork(3) = vhaes_aux%get_ldwork(nlat)
        ldwork(4) = vhses_aux%get_ldwork(nlat)

        return_value = maxval(ldwork)

    end function get_ldwork

    pure subroutine get_legendre_workspace(nlat, nlon, workspace, nt, ityp)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: workspace(:)
        integer(ip), optional, intent(in)  :: nt
        integer(ip), optional, intent(in)  :: ityp

        ! Local variables
        type(ShaesAux) :: shaes_aux
        type(ShsesAux) :: shses_aux
        type(VhaesAux) :: vhaes_aux
        type(VhsesAux) :: vhses_aux
        integer(ip)    :: work_size(4)
        integer(ip)    :: lwork, nt_op, ityp_op

        !  Address optional arguments
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

        ! Get required workspace size
        work_size(1) = shaes_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(2) = shses_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(3) = vhaes_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        work_size(4) = vhses_aux%get_legendre_workspace_size(nlat, nlon, nt_op, ityp_op)
        lwork = maxval(work_size)

        !  Allocate memory
        allocate( workspace(lwork) )

    end subroutine get_legendre_workspace

end module type_RegularWorkspace
