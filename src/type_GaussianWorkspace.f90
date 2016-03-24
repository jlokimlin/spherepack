module type_GaussianWorkspace

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_Workspace, only: &
        Workspace

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: GaussianWorkspace


    ! Declare derived data type
    type, extends (Workspace), public :: GaussianWorkspace
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure,         public  :: create => create_gaussian_workspace
        procedure,         public  :: destroy => destroy_gaussian_workspace
        procedure,         private :: initialize_gaussian_scalar_analysis
        procedure,         private :: initialize_gaussian_scalar_synthesis
        procedure,         private :: initialize_gaussian_scalar_transform
        procedure,         private :: initialize_gaussian_vector_analysis
        procedure,         private :: initialize_gaussian_vector_synthesis
        procedure,         private :: initialize_gaussian_vector_transform
        procedure, nopass, private :: get_lshags
        procedure, nopass, private :: get_lshsgs
        procedure, nopass, private :: get_lvhags
        procedure, nopass, private :: get_lvhsgs
        final                      :: finalize_gaussian_workspace
        !----------------------------------------------------------------------
    end type GaussianWorkspace


contains


    subroutine create_gaussian_workspace( this, nlat, nlon )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        ! Set up scalar transform
        call this%initialize_gaussian_scalar_transform(nlat, nlon)

        ! Set up vector transform
        call this%initialize_gaussian_vector_transform(nlat, nlon)

        ! Set flag
        this%initialized = .true.

    end subroutine create_gaussian_workspace

    subroutine destroy_gaussian_workspace( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if ( this%initialized .eqv. .false.) return

        ! Release memory from parent type
        call this%destroy_workspace()

        ! Set flag
        this%initialized = .true.

    end subroutine destroy_gaussian_workspace


    subroutine initialize_gaussian_scalar_analysis( this, nlat, nlon )
        !
        ! Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lshags
        real (wp), allocatable :: dwork(:)
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        lwork = this%get_lwork(nlat, nlon)
        ldwork = this%get_ldwork(nlat)
        lshags = this%get_lshags(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%legendre_workspace)) deallocate( this%legendre_workspace )
        if (allocated(this%forward_scalar)) deallocate( this%forward_scalar )

        ! Allocate memory
        allocate( this%legendre_workspace(lwork) )
        allocate( this%forward_scalar(lshags) )
        allocate( dwork(lwork) )

        ! Compute workspace arrays
        associate( &
            wshags => this%forward_scalar, &
            work => this%legendre_workspace, &
            ierror => error_flag &
            )
            call shagsi( nlat, nlon, wshags, lshags, work, lwork, dwork, ldwork, ierror )
        end associate

        ! Release memory
        deallocate( dwork )


        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE (GaussianWorkspace): '&
                    //'Error in the specification of NLAT'
            case(2)
                error stop 'TYPE (GaussianWorkspace): '&
                    //'Error in the specification of NLON'
            case(3)
                error stop 'TYPE (GaussianWorkspace): '&
                    //'Error in the specification of extent for FORWARD_SCALAR'
            case(4)
                error stop 'TYPE (GaussianWorkspace): '&
                    //'Error in the specification of extent for LEGENDRE_WORKSPACE'
            case(5)
                error stop 'TYPE (GaussianWorkspace): '&
                    //'Error in the specification of extent for DWORK'
            case(6)
                error stop 'TYPE (GaussianWorkspace): '&
                    //'Error in call to GAQD to compute gaussian points: failure in eigenvalue routine)'
            case default
                error stop 'TYPE (GaussianWorkspace): '&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_scalar_analysis



    subroutine initialize_gaussian_scalar_synthesis( this, nlat, nlon )
        !
        ! Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic synthesis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#shags.html
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: lwork, ldwork, lshsgs
        real (wp), allocatable :: dwork(:)
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        lwork = this%get_lwork(nlat, nlon)
        ldwork = this%get_ldwork(nlat)
        lshsgs = this%get_lshsgs(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%legendre_workspace)) deallocate( this%legendre_workspace )
        if (allocated(this%backward_scalar)) deallocate( this%backward_scalar )

        ! Allocate memory
        allocate( this%legendre_workspace(lwork) )
        allocate( dwork(ldwork) )
        allocate( this%backward_scalar(lshsgs) )

        ! Compute workspace
        associate( &
            wshsgs => this%backward_scalar, &
            work => this%legendre_workspace, &
            ierror => error_flag &
            )
            call shsgsi(nlat, nlon, wshsgs, lshsgs, work, lwork, dwork, ldwork, ierror)
        end associate

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_SCALAR_SYNTHESIS'&
                    //'Error in the specification of NLAT'
            case(2)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_SCALAR_SYNTHESIS'&
                    //'Error in the specification of NLON'
            case(3)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_SCALAR_SYNTHESIS'&
                    //'Error in the specification of extent for BACKWARD_SCALAR'
            case(4)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_SCALAR_SYNTHESIS'&
                    //'Error in the specification of extent for LEGENDRE_WORKSPACE'
            case(5)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_SCALAR_SYNTHESIS'&
                    //'Error in the specification of extent for DWORK'
            case(6)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_SCALAR_SYNTHESIS'&
                    //'Error in call to GAQD: due to failure in eigenvalue routine'
            case default
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_SCALAR_SYNTHESIS'&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_scalar_synthesis


    subroutine initialize_gaussian_scalar_transform( this, nlat, nlon )
        !
        ! Purpose:
        !
        !  Set the various workspace arrays and pointers to perform the
        !  (real) scalar harmonic transform on a gaussian grid
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Set up scalar analysis
        call this%initialize_gaussian_scalar_analysis(nlat, nlon)

        ! Set up scalar synthesis
        call this%initialize_gaussian_scalar_synthesis(nlat, nlon)

        ! Allocate memory
        allocate( this%real_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_harmonic_coefficients(nlat, nlat) )

    end subroutine initialize_gaussian_scalar_transform


    subroutine initialize_gaussian_vector_analysis( this, nlat, nlon )
        !
        ! Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#vhags.html
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: ldwork, lvhags
        real (wp), allocatable :: dwork(:)
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        ldwork = this%get_ldwork(nlat)
        lvhags = this%get_lvhags(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%forward_vector)) deallocate( this%forward_vector )

        ! Allocate memory
        allocate( dwork(ldwork) )
        allocate( this%forward_vector(lvhags) )

        ! Compute workspace
        associate( &
            wvhags => this%forward_vector, &
            ierror => error_flag &
            )
            call vhagsi( nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror )
        end associate

        ! Release memory
        deallocate( dwork )

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of NLAT'
            case(2)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of NLON'
            case(3)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of extent for FORWARD_VECTOR'
            case(4)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Error in the specification of extent for DWORK'
            case default
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_GAUSSIAN_VECTOR_ANALYSIS'&
                    //'Undetermined error flag'
        end select

    end subroutine initialize_gaussian_vector_analysis


    subroutine initialize_gaussian_vector_synthesis( this, nlat, nlon )
        !
        !< Purpose:
        !
        !  Set the various workspace arrays to perform the
        !  (real) scalar harmonic analysis on a gaussian grid
        !
        !  Reference:
        !  https://www2.cisl.ucar.edu/spherepack/documentation#vhsgsi.html
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)           :: error_flag
        integer (ip)           :: ldwork, lvhsgs
        real (wp), allocatable :: dwork(:)
        !----------------------------------------------------------------------

        ! Compute dimensions of various workspace arrays
        ldwork = this%get_ldwork(nlat)
        lvhsgs = this%get_lvhsgs(nlat, nlon)

        ! Release memory ( if necessary )
        if (allocated(this%backward_vector) ) deallocate( this%backward_vector )

         ! Allocate memory
        allocate( dwork(ldwork) )
        allocate( this%backward_vector(lvhsgs) )

        ! Compute workspace
        associate( &
            wvhsgs => this%backward_vector, &
            ierror => error_flag &
            )
            call vhsgsi( nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror )
        end associate

        ! Release memory
        deallocate( dwork )

        ! Address error flag
        select case (error_flag)
            case(0)
                return
            case(1)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_VECTOR_synthesis'&
                    //'Error in the specification of NLAT'
            case(2)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_VECTOR_SYNTHESIS'&
                    //'Error in the specification of NLON'
            case(3)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_VECTOR_SYNTHESIS'&
                    //'Error in the specification of extent for BACKWARD_VECTOR'
            case(4)
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_VECTOR_SYNTHESIS'&
                    //'Error in the specification of extent for DWORK'
            case default
                error stop 'TYPE (GaussianWorkspace) in INITIALIZE_VECTOR_SYNTHESIS'&
                    //'Undetermined error flag'
        end select


    end subroutine initialize_gaussian_vector_synthesis


    subroutine initialize_gaussian_vector_transform( this, nlat, nlon )
        !
        ! Purpose:
        !
        !  Sets the various workspace arrays and pointers
        !  required for the (real) vector harmonic transform
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (GaussianWorkspace), intent (in out) :: this
        integer (ip),              intent (in)     :: nlat
        integer (ip),              intent (in)     :: nlon
        !----------------------------------------------------------------------

        ! Set up vector analysis
        call this%initialize_gaussian_vector_analysis(nlat, nlon)

        ! Set up vector analysis
        call this%initialize_gaussian_vector_synthesis(nlat, nlon)

        ! Allocate memory
        allocate( this%real_polar_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_polar_harmonic_coefficients(nlat, nlat) )
        allocate( this%real_azimuthal_harmonic_coefficients(nlat, nlat) )
        allocate( this%imaginary_azimuthal_harmonic_coefficients(nlat, nlat) )

    end subroutine initialize_gaussian_vector_transform


    pure function get_lshags(nlat, nlon) result ( return_value )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)::  l1, l2
        !----------------------------------------------------------------------

        ! Compute parity
        if ( mod(nlon,2) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod(nlat,2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            nlat * (3 * (l1 + l2) - 2) &
            +(l1 - 1) * (l2 * (2 * nlat - l1) - 3 * l1)/2 &
            + nlon + 15

    end function get_lshags


    pure function get_lshsgs(nlat, nlon) result ( return_value )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        integer (ip)              :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: l1,  l2
        !----------------------------------------------------------------------

        ! Compute parity
        if ( mod(nlon,2) == 0 ) then
            l1 = min( nlat, (nlon + 2)/2 )
        else
            l1 = min( nlat, (nlon + 1)/2 )
        end if

        if ( mod(nlat,2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value = &
            nlat * (3 * (l1 + l2) - 2) &
            +(l1 - 1) * (l2 * (2 * nlat - l1) - 3 * l1)/2 &
            + nlon + 15

    end function get_lshsgs


    pure function get_lvhags(nlat, nlon) result ( return_value )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        integer (ip)              :: return_value
        !----------------------------------------------------------------------

        return_value = &
            (nlat + 1) * (nlat + 1) * nlat / 2 &
            + nlon + 15

    end function get_lvhags


    pure function get_lvhsgs(nlat, nlon) result ( return_value )
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
        if ( mod(nlon,2) == 0 ) then
            l1 = min( nlat, nlon/2 )
        else
            l1 = min( nlat,(nlon + 1)/2 )
        end if

        if ( mod(nlat,2) == 0 ) then
            l2 = nlat/2
        else
            l2 = (nlat + 1)/2
        end if

        return_value =  &
            l1 * l2 * (nlat + nlat - l1 + 1) &
            + nlon + 15 + 2 * nlat

    end function get_lvhsgs


    subroutine finalize_gaussian_workspace( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        type (GaussianWorkspace), intent (in out) :: this
        !----------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_gaussian_workspace

end module type_GaussianWorkspace
