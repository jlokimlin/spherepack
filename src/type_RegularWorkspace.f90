module type_RegularWorkspace

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_Workspace, only: &
        Workspace

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
        final                      :: finalize_regular_workspace
        !----------------------------------------------------------------------
    end type RegularWorkspace


contains


    subroutine create_regular_workspace( this, nlat, nlon )
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



    subroutine destroy_regular_workspace( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (RegularWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if ( this%initialized .eqv. .false. ) return

        ! Release memory from parent type
        call this%destroy_workspace()

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_regular_workspace


    subroutine initialize_regular_scalar_analysis( this, nlat, nlon )
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
        lwork = this%get_lwork(nlat, nlon)
        ldwork = this%get_ldwork(nlat)
        lshaes = this%get_lshaes(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%legendre_workspace)) deallocate( this%legendre_workspace )
        if (allocated(this%forward_scalar)) deallocate( this%forward_scalar )

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
            call shaesi( nlat, nlon, wshaes, lshaes, work, lwork, dwork, ldwork, ierror )
        end associate

        ! Release memory
        deallocate( dwork )

         !  Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_ANALYSIS'&
                    //'Error in the specification of extent for FORWARD_SCALAR'
            case(4)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_ANALYSIS'&
                    //'Error in the specification of extent for LEGENDRE_WORKSPACE'
            case default
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_ANALYSIS'&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_regular_scalar_analysis



    subroutine initialize_regular_scalar_synthesis( this, nlat, nlon )
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
        lwork = this%get_lwork(nlat, nlon)
        ldwork = this%get_ldwork(nlat)
        lshses = this%get_lshses(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%legendre_workspace)) deallocate( this%legendre_workspace )
        if (allocated(this%backward_scalar)) deallocate( this%backward_scalar )

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
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_SYNTHESIS'&
                    //'Error in the specification of extent for FORWARD_SCALAR'
            case(4)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_SYNTHESIS'&
                    //'Error in the specification of extent for LEGENDRE_WORKSPACE'
            case default
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_SCALAR_SYNTHESIS'&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_regular_scalar_synthesis



    subroutine initialize_regular_scalar_transform( this, nlat, nlon )
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



    subroutine initialize_regular_vector_analysis( this, nlat, nlon )
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
        if (allocated(this%forward_vector)) deallocate( this%forward_vector )

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
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_ANALYSIS'&
                    //'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_ANALYSIS'&
                    //' Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_ANALYSIS'&
                    //' Error in the specification of extent for FORWARD_VECTOR'
            case(4)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_ANALYSIS'&
                    //' Error in the specification of extent for unsaved WORK'
            case(5)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_ANALYSIS'&
                    //' Error in the specification of extent for unsaved DWORK'
            case default
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_ANALYSIS'&
                    //' Undetermined error flag'
        end select


    end subroutine initialize_regular_vector_analysis



    subroutine initialize_regular_vector_synthesis( this, nlat, nlon )
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
        if (allocated(this%backward_vector)) deallocate( this%backward_vector )

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
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NUMBER_OF_LATITUDES'
            case(2)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_SYNTHESIS'&
                    //' Error in the specification of NUMBER_OF_LONGITUDES'
            case(3)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_SYNTHESIS'&
                    //' Error in the specification of extent for BACKWARD_VECTOR'
            case(4)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_SYNTHESIS'&
                    //' Error in the specification of extent for unsaved WORK'
            case(5)
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_SYNTHESIS'&
                    //' Error in the specification of extent for unsaved DWORK'
            case default
                error stop 'TYPE (RegularWorkspace) in INITIALIZE_REGULAR_VECTOR_SYNTHESIS'&
                    //' Undetermined error flag'
        end select


    end subroutine initialize_regular_vector_synthesis



    subroutine initialize_regular_vector_transform( this, nlat, nlon )
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
        allocate( this%imaginary_azimuthal_harmonic_coefficients( nlat, nlat) )

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
        if ( mod(nlon, 2) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod(nlat, 2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            (l1 * l2 * (nlat + nlat - l1 + 1))/2 + nlon + 15

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
        if ( mod(nlon, 2) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod(nlat, 2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            (l1 * l2 * (nlat + nlat - l1 + 1))/2 + nlon + 15

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
        if ( mod(nlon, 2) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod(nlat, 2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            l1 * l2 * (nlat + nlat - l1 + 1) + nlon + 15

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
        if ( mod(nlon, 2) == 0 ) then
            l1 = min( nlat, nlon/2 )
        else
            l1 = min( nlat,(nlon + 1)/2 )
        end if

        if ( mod(nlat, 2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value =  &
            l1 * l2 * (nlat + nlat - l1 + 1) + nlon + 15

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
        if ( mod(nlon, 2) == 0 ) then
            l1 = min( nlat, nlon/2 )
        else
            l1 = min( nlat,(nlon + 1)/2 )
        end if

        if ( mod(nlat, 2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = 3 * (max(l1-2,0) * (nlat+nlat-l1-1))/2 + 5*l2*nlat

    end function get_lwork_unsaved


    subroutine finalize_regular_workspace( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (RegularWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy_workspace()

    end subroutine finalize_regular_workspace


end module type_RegularWorkspace
