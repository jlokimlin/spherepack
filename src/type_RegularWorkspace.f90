module type_RegularWorkspace

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_Workspace, only: &
        Workspace

    use module_shaes, only: &
        shaesi

    use module_shses, only: &
        shsesi

    use module_vhaes, only: &
        vhaesi

    use module_vhses, only: &
        vhsesi

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
        procedure,         public  :: create => create_regular_workspace
        procedure,         public  :: destroy => destroy_regular_workspace
        procedure,         private :: initialize_regular_scalar_analysis
        procedure,         private :: initialize_regular_scalar_synthesis
        procedure,         private :: initialize_regular_scalar_transform
        procedure,         private :: initialize_regular_vector_analysis
        procedure,         private :: initialize_regular_vector_synthesis
        procedure,         private :: initialize_regular_vector_transform
        procedure, nopass, private :: get_lshaes
        procedure, nopass, private :: get_lshses
        procedure, nopass, private :: get_lvhaes
        procedure, nopass, private :: get_lvhses
        procedure, nopass, private :: get_lwork_unsaved
        generic,           public  :: assignment (=) => copy_regular_workspace
        procedure,         private :: copy_regular_workspace
        final                      :: finalize_regular_workspace
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
        integer (ip),         intent (in) :: nlat !! number of latitudinal points 0 <= theta <= pi
        integer (ip),         intent (in) :: nlon !! number of longitudinal points 0 <= phi <= 2*pi
        type (RegularWorkspace)           :: return_value
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
        if (this%initialized .eqv. .false.) then
            return
        end if

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
        real (wp), allocatable :: dwork(:)
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        lwork = max(this%get_lwork(nlat, nlon), 5*(nlat**2)*nlon)
        ldwork = max(this%get_ldwork(nlat), 12*nlat*nlon, 4*(nlat**2))
        lshaes = this%get_lshaes(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%legendre_workspace)) then
            deallocate( this%legendre_workspace )
        end if

        if (allocated(this%forward_scalar)) then
            deallocate( this%forward_scalar )
        end if

        ! Allocate memory
        allocate( this%legendre_workspace(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%forward_scalar(lshaes) )

        ! Compute workspace
        associate( &
            wshaes => this%forward_scalar, &
            work => this%legendre_workspace, &
            ierror => error_flag &
            )
            call shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, ldwork, ierror)
        end associate

        ! Release memory
        deallocate( dwork )

         !  Address error flag
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

    end subroutine initialize_regular_scalar_analysis



    subroutine initialize_regular_scalar_synthesis(this, nlat, nlon)
        !
        !  Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a regular (equally-spaced) grid.
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shsesi.html
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
        integer (ip)           :: lwork, ldwork, lshses
        real (wp), allocatable :: dwork(:)
        !----------------------------------------------------------------------

        ! Set up various workspace dimensions
        lwork = max(this%get_lwork(nlat, nlon), 5*(nlat**2)*nlon)
        ldwork = max(this%get_ldwork(nlat), 12*nlat*nlon, 4*(nlat**2))
        lshses = this%get_lshses(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%legendre_workspace)) then
            deallocate( this%legendre_workspace )
        end if

        if (allocated(this%backward_scalar)) then
            deallocate( this%backward_scalar )
        end if

        ! Allocate memory
        allocate( this%legendre_workspace(lwork) )
        allocate( this%backward_scalar(lshses) )
        allocate( dwork(ldwork) )

        ! Compute workspace
        associate( &
            wshses => this%backward_scalar, &
            work => this%legendre_workspace, &
            ierror => error_flag &
            )
            call shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, ldwork, ierror)
        end associate

        ! Release memory
        deallocate( dwork )

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
        !
        !< Purpose:
        !
        ! Sets the various workspace arrays required for
        ! (real) vector harmonic analysis on a regular (equally-spaced) grid
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhaesi.html
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
        integer (ip)           :: lwork, ldwork, lvhaes
        real (wp), allocatable :: work(:), dwork(:)
        !----------------------------------------------------------------------

        ! Compute various workspace dimensions
        lwork = this%get_lwork_unsaved(nlat, nlon)
        ldwork = this%get_ldwork(nlat)
        lvhaes = this%get_lvhaes(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%forward_vector)) then
            deallocate( this%forward_vector )
        end if

        ! Allocate memory
        allocate( work(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%forward_vector(lvhaes) )

        ! Compute workspace
        associate( &
            wvhaes => this%forward_vector, &
            ierror => error_flag &
            )
            call vhaesi(nlat, nlon, wvhaes, lvhaes, work, lwork, dwork, ldwork, ierror)
        end associate

        ! Release memory
        deallocate( work )
        deallocate( dwork )

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


    end subroutine initialize_regular_vector_analysis



    subroutine initialize_regular_vector_synthesis(this, nlat, nlon)
        !
        ! Purpose:
        !
        ! Sets the various workspace arrays required for
        ! (real) vector harmonic analysis on a regular (equally-spaced) grid
        !
        ! Reference:
        ! https://www2.cisl.ucar.edu/spherepack/documentation#vhsesi.html
        !
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
        integer (ip)           :: lwork, ldwork, lvhses
        real (wp), allocatable :: work(:), dwork(:)
        !----------------------------------------------------------------------

        ! Compute various workspace dimensions
        lwork = this%get_lwork_unsaved(nlat, nlon)
        ldwork = this%get_ldwork(nlat)
        lvhses = this%get_lvhses(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%backward_vector)) then
            deallocate( this%backward_vector )
        end if

        ! Allocate memory
        allocate( work(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%backward_vector(lvhses) )

        ! Compute workspace
        associate( &
            wvhses => this%backward_vector, &
            ierror => error_flag &
            )
            call vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, ldwork, ierror)
        end associate

        ! Release memory
        deallocate( work )
        deallocate( dwork )

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



    pure function get_lshaes(nlat, nlon) result ( return_value )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: l1,  l2
        !----------------------------------------------------------------------

        ! Compute parity
        if (mod(nlon, 2) == 0) then
            l1 = min(nlat, (nlon+2)/2)
        else
            l1 = min(nlat, (nlon+1)/2)
        end if

        if (mod(nlat, 2) == 0) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = ( l1 * l2 * (2*nlat-l1+1) )/2 + nlon+15

    end function get_lshaes



    pure function get_lshses(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        integer (ip)              :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: l1, l2
        !----------------------------------------------------------------------

        ! Compute parity
        if (mod(nlon, 2) == 0) then
            l1 = min(nlat, (nlon+2)/2)
        else
            l1 = min(nlat, (nlon+1)/2)
        end if

        if (mod(nlat, 2) == 0) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = get_lshaes(nlat, nlon)


    end function get_lshses


    pure function get_lvhaes(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        integer (ip) :: l1, l2
        !----------------------------------------------------------------------

        ! Compute parity
        if (mod(nlon, 2) == 0) then
            l1 = min(nlat, (nlon+2)/2)
        else
            l1 = min(nlat, (nlon+1)/2)
        end if

        if (mod(nlat, 2) == 0) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = l1 * l2 * (2*nlat-l1+1) + nlon+15

    end function get_lvhaes



    pure function get_lvhses(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: l1, l2
        !----------------------------------------------------------------------

        ! Compute parity
        if (mod(nlon, 2) == 0) then
            l1 = min(nlat, nlon/2 )
        else
            l1 = min(nlat, (nlon+1)/2)
        end if

        if (mod(nlat, 2) == 0) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = l1 * l2 * (2*nlat-l1+1) + nlon+15

    end function get_lvhses



    pure function get_lwork_unsaved(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: l1, l2
        !----------------------------------------------------------------------

        ! Compute parity
        if (mod(nlon, 2) == 0) then
            l1 = min(nlat, nlon/2 )
        else
            l1 = min(nlat, (nlon+1)/2)
        end if

        if (mod(nlat, 2) == 0) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        associate( &
            a => 3 * (max(l1-2, 0) * (nlat+nlat-l1-1))/2 + 5*l2*nlat, &
            b => 5 * nlat * l2 + 3 * ( (l1-2)*(nlat+nlat-l1-1) )/2 &
            )
            return_value = max(a, b)
        end associate

    end function get_lwork_unsaved



    subroutine finalize_regular_workspace(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (RegularWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy_workspace()

    end subroutine finalize_regular_workspace



end module type_RegularWorkspace
