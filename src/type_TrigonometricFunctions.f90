
!
!< Author:
! Jon Lo Kim Lin
!

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

    !----------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------------------------------------
    character (len=250) :: error_message     !! Probably long enough
    integer (ip)        :: allocate_status   !! To check allocation status
    integer (ip)        :: deallocate_status !! To check deallocation status
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, public :: TrigonometricFunctions
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                 public :: initialized = .false. !! Initialization flag
        integer (ip),            public :: NLON = 0 !! number of longitudinal points in phi
        integer (ip),            public :: NLAT = 0 !! number of latitudinal points in theta
        real (wp), allocatable, public :: sint(:)  !! sin(theta): 0 <= theta <= pi
        real (wp), allocatable, public :: cost(:)  !! cos(theta): 0 <= theta <= pi
        real (wp), allocatable, public :: sinp(:)  !! sin(phi):   0 <=  phi  <= 2*pi
        real (wp), allocatable, public :: cosp(:)  !! cos(phi):   0 <=  phi  <= 2*pi
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public :: create => create_trigonometric_functions
        procedure, public :: destroy => destroy_trigonometric_functions
        final              :: finalize_trigonometric_functions
        !----------------------------------------------------------------------
    end type TrigonometricFunctions

contains
    !
    
    !
    subroutine create_trigonometric_functions( this, latitudinal_grid, longitudinal_grid )
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

            allocate( &
                this%sint( nlat ), &
                this%cost( nlat ), &
                this%sinp( nlon ), &
                this%cosp( nlon ), &
                stat=allocate_status, &
                errmsg = error_message )

            ! Check allocation status
            if ( allocate_status /= 0 ) then
                write( stderr, '(A)') 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)') 'Allocation failed in create_trigonometric_functions'
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

    end subroutine create_trigonometric_functions
    !
    
    !
    subroutine destroy_trigonometric_functions( this )
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
        ! Release memory
        !--------------------------------------------------------------------------------

         ! Check if array is allocated
        if ( allocated( this%sint ) ) then

            ! Deallocate array
            deallocate( &
                this%sint, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating SINT failed in destroy_trigonometric_functions'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%cost ) ) then

            ! Deallocate array
            deallocate( &
                this%cost, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating COST failed in destroy_trigonometric_functions'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%sinp ) ) then

            ! Deallocate array
            deallocate( &
                this%sinp, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating SINP failed in destroy_trigonometric_functions'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%cosp ) ) then

            ! Deallocate array
            deallocate( &
                this%cosp, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (TrigonometricFunctions)'
                write( stderr, '(A)' ) 'Deallocating COSP failed in destroy_trigonometric_functions'
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

    end subroutine destroy_trigonometric_functions
    !
    
    !
    subroutine finalize_trigonometric_functions( this )
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (TrigonometricFunctions), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_trigonometric_functions
    !
    
    !
end module type_TrigonometricFunctions
