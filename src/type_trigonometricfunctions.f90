!*****************************************************************************************
!
!< Author:
! Jon Lo Kim Lin
!
!*****************************************************************************************
!
module type_TrigonometricFunctions

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: TrigonometricFunctions

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    character (len=250) :: error_message     !! Probably long enough
    integer (ip)        :: allocate_status   !! To check allocation status
    integer (ip)        :: deallocate_status !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, public :: TrigonometricFunctions

        ! All components are public unless stated otherwise

        logical                           :: initialized = .false. !! Instantiation status
        integer (ip)                      :: NLON        = 0        !! number of longitudinal points
        integer (ip)                      :: NLAT        = 0        !! number of latitudinal points
        real (wp), allocatable           :: sint(:)                !! sin(theta): 0 <= theta <= pi
        real (wp), allocatable           :: cost(:)                !! cos(theta): 0 <= theta <= pi
        real (wp), allocatable           :: sinp(:)                !! sin(phi):   0 <=  phi  <= 2*pi
        real (wp), allocatable           :: cosp(:)                !! cos(phi):   0 <=  phi  <= 2*pi

    contains

        ! All method are private unless stated otherwise
        private

        !---------------------------------------------------------------------------------
        ! Public methods
        !---------------------------------------------------------------------------------
        procedure,                    public   :: create => create_trigonometricfunctions
        procedure,                    public   :: destroy => destroy_trigonometricfunctions
        !---------------------------------------------------------------------------------
        ! Finalizer
        !---------------------------------------------------------------------------------
        final                                :: finalize_trigonometricfunctions
        !---------------------------------------------------------------------------------

    end type TrigonometricFunctions

contains
    !
    !*****************************************************************************************
    !
    subroutine create_trigonometricfunctions( this, latitudinal_grid, longitudinal_grid )
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (TrigonometricFunctions), intent (in out) :: this
        real (wp),                      intent (in)     :: latitudinal_grid(:)
        real (wp),                      intent (in)     :: longitudinal_grid(:)
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized ) then
            call this%destroy()
        end if

        !--------------------------------------------------------------------------------
        ! Set constants
        !--------------------------------------------------------------------------------

        this%NLAT = size( latitudinal_grid )
        this%NLON = size( longitudinal_grid )

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            allocate ( &
                this%sint( nlat ), &
                this%cost( nlat ), &
                this%sinp( nlon ), &
                this%cosp( nlon ), &
                stat   = allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then
                write( stderr, '(A)') 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)') 'Allocation failed in CREATE_TRIGONOMETRICFUNCTIONS'
                write( stderr, '(A)') trim( error_message )
            end if

        end associate

        !--------------------------------------------------------------------------------
        ! compute trigonometric functions
        !--------------------------------------------------------------------------------

        associate( &
            theta => latitudinal_grid, &
            phi   => longitudinal_grid, &
            sint  => this%sint, &
            cost  => this%cost, &
            sinp  => this%sinp, &
            cosp  => this%cosp &
            )

            sint = sin( theta )
            cost = cos( theta )
            sinp = sin( phi )
            cosp = cos( phi )

        end associate

        !--------------------------------------------------------------------------------
        ! Set initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .true.

    end subroutine create_trigonometricfunctions
    !
    !*****************************************************************************************
    !
    subroutine destroy_trigonometricfunctions( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (TrigonometricFunctions), intent (in out) :: this
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check status
        !--------------------------------------------------------------------------------

        if ( .not. this%initialized ) return

        !--------------------------------------------------------------------------------
        ! Deallocate arrays
        !--------------------------------------------------------------------------------

         ! Check if array is allocated
        if ( allocated( this%sint ) ) then

            ! Deallocate array
            deallocate( &
                this%sint, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating SINT failed in DESTROY_TRIGONOMETRICFUNCTIONS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%cost ) ) then

            ! Deallocate array
            deallocate( &
                this%cost, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating COST failed in DESTROY_TRIGONOMETRICFUNCTIONS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%sinp ) ) then

            ! Deallocate array
            deallocate( &
                this%sinp, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating SINP failed in DESTROY_TRIGONOMETRICFUNCTIONS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%cosp ) ) then

            ! Deallocate array
            deallocate( &
                this%cosp, &
                stat   = deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating COSP failed in DESTROY_TRIGONOMETRICFUNCTIONS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        !--------------------------------------------------------------------------------
        ! Reset constants
        !--------------------------------------------------------------------------------

        this%NLON  = 0
        this%NLAT  = 0

        !--------------------------------------------------------------------------------
        ! Reset initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .false.

    end subroutine destroy_trigonometricfunctions
    !
    !*****************************************************************************************
    !
    subroutine finalize_trigonometricfunctions( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (TrigonometricFunctions), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_trigonometricfunctions
    !
    !*****************************************************************************************
    !
end module type_TrigonometricFunctions
