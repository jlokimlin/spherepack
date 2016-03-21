
!
!< Author:
! Jon Lo Kim Lin
!

!
module type_SphericalUnitVectors

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stderr => ERROR_UNIT

    use type_ThreeDimensionalVector, only: &
        ThreeDimensionalVector, &
        assignment(=), &
        operator(*)

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: SphericalUnitVectors

    !----------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------------------------------------
    character (len=250) :: error_message     !! Probably long enough
    integer (ip)        :: allocate_status   !! To check allocation status
    integer (ip)        :: deallocate_status !! To check deallocation status
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, public :: SphericalUnitVectors

        ! All components are public unless stated otherwise

        logical                                      :: initialized = .false. !! Instantiation status
        integer (ip)                                 :: NLON        = 0        !! number of longitudinal points
        integer (ip)                                 :: NLAT        = 0        !! number of latitudinal points
        type (ThreeDimensionalVector), allocatable  :: radial(:, :)
        type (ThreeDimensionalVector), allocatable  :: polar(:, :)
        type (ThreeDimensionalVector), allocatable  :: azimuthal(:, :)

    contains

        ! All method are private unless stated otherwise
        private

        !----------------------------------------------------------------------
        ! Methods
        !----------------------------------------------------------------------
        procedure, public   :: create => create_sphericalunitvectors
        procedure, public   :: destroy => destroy_sphericalunitvectors
        procedure, public   :: get_spherical_angle_components
        procedure            :: assert_initialized
        final                :: finalize_sphericalunitvectors
        !----------------------------------------------------------------------

    end type SphericalUnitVectors

contains
    !
    
    !
    subroutine create_sphericalunitvectors( this, sint, cost, sinp, cosp )
        !
        !< Purpose:
        ! Sets the spherical unit vectors
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out) :: this
        real (wp),                    intent (in)     :: sint(:)
        real (wp),                    intent (in)     :: cost(:)
        real (wp),                    intent (in)     :: sinp(:)
        real (wp),                    intent (in)     :: cosp(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)            ::  k,  l !! Counters
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized ) then
            call this%destroy()
        end if

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        this%NLAT = size( sint )
        this%NLON = size( cosp )

        !--------------------------------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON  &
            )

            ! Allocate arrays
            allocate( &
                this%radial(    nlat, nlon ), &
                this%polar(     nlat, nlon ), &
                this%azimuthal( nlat, nlon ), &
                stat=allocate_status, &
                errmsg = error_message )

            ! Check allocate status
            if ( allocate_status /= 0 ) then
                write( stderr, '(A)') 'TYPE (SphericalUnitVectors)'
                write( stderr, '(A)') 'Allocation failed in CREATE_SPHERICALUNITVECTORS'
                write( stderr, '(A)') trim( error_message )
            end if
        end associate

        !--------------------------------------------------------------------------------
        ! Compute spherical unit vectors
        !--------------------------------------------------------------------------------

        associate( &
            nlat  => this%NLAT, &
            nlon  => this%NLON, &
            r     => this%radial, &
            theta => this%polar, &
            phi   => this%azimuthal &
            )

            ! compute spherical unit vectors
            do l = 1, nlon
                do k = 1, nlat

                    ! set radial unit vector
                    r(k, l) = &
                        ThreeDimensionalVector( &
                        x = sint(k) * cosp(l), &
                        y = sint(k) * sinp(l), &
                        z = cost(k) &
                        )

                    ! set polar unit vector
                    theta(k, l) = &
                        ThreeDimensionalVector( &
                        x = cost(k) * cosp(l), &
                        y = cost(k) * sinp(l), &
                        z = -sint(k) &
                        )

                    ! set azimuthal unit vector
                    phi(k, l) = &
                        ThreeDimensionalVector( &
                        x = -sinp(l), &
                        y =  cosp(l), &
                        z =  0.0_wp &
                        )
                end do
            end do
        end associate

        !--------------------------------------------------------------------------------
        ! Set initialization flag
        !--------------------------------------------------------------------------------

        this%initialized = .true.

    end subroutine create_sphericalunitvectors
    !
    
    !
    subroutine destroy_sphericalunitvectors( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out) :: this
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check status
        !--------------------------------------------------------------------------------

        if ( .not. this%initialized ) return

        !--------------------------------------------------------------------------------
        ! Release memory
        !--------------------------------------------------------------------------------

        ! Check if array is allocated
        if ( allocated( this%radial ) ) then

            ! Deallocate array
            deallocate( &
                this%radial, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (SphericalUnitVectors)'
                write( stderr, '(A)' ) 'Deallocating RADIAL failed in DESTROY_SPHERICALUNITVECTORS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%polar ) ) then

            ! Deallocate array
            deallocate( &
                this%polar, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (SphericalUnitVectors)'
                write( stderr, '(A)' ) 'Deallocating POLAR failed in DESTROY_SPHERICALUNITVECTORS'
                write( stderr, '(A)' ) trim( error_message )
            end if
        end if

        ! Check if array is allocated
        if ( allocated( this%azimuthal ) ) then

            ! Deallocate array
            deallocate( &
                this%azimuthal, &
                stat=deallocate_status, &
                errmsg = error_message )

            ! Check deallocation status
            if ( deallocate_status /= 0 ) then
                write( stderr, '(A)' ) 'TYPE (SphericalUnitVectors)'
                write( stderr, '(A)' ) 'Deallocating AZIMUTHAL failed in DESTROY_SPHERICALUNITVECTORS'
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

    end subroutine destroy_sphericalunitvectors
    !
    
    !
    subroutine assert_initialized( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. this%initialized ) then

            write( stderr, '(A)' ) 'TYPE (SphericalUnitVectors)'
            write( stderr, '(A)' ) 'You must instantiate object before calling methods'

        end if

    end subroutine assert_initialized
    !
    
    !
    subroutine get_spherical_angle_components( this, &
        vector_function, polar_component, azimuthal_component )
        !
        !< Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (SphericalUnitVectors), intent (in out)  :: this
        real (wp),                    intent (in)      :: vector_function(:, :, :)
        real (wp),                    intent (out)     :: polar_component(:, :)
        real (wp),                    intent (out)     :: azimuthal_component(:, :)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip)                  :: k, l         !! Counters
        type (ThreeDimensionalVector) :: vector_field !! To convert array to vector
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Check if object is usable
        !--------------------------------------------------------------------------------

        call this%assert_initialized()

        !--------------------------------------------------------------------------------
        ! Calculate the spherical angle components
        !--------------------------------------------------------------------------------

        associate( &
            nlat => this%NLAT, &
            nlon => this%NLON &
            )

            do l = 1, nlon
                do k = 1, nlat

                    ! Convert array to vector
                    vector_field = vector_function(:, k, l)

                    associate( &
                        theta => this%polar( k, l ), &
                        phi   => this%azimuthal( k, l ) &
                        )

                        ! set the theta component
                        polar_component( k, l ) = theta.dot.vector_field

                        ! set the azimuthal_component
                        azimuthal_component( k, l ) = phi.dot.vector_field

                    end associate
                end do
            end do
        end associate

    end subroutine get_spherical_angle_components
    !
    
    !
    subroutine finalize_sphericalunitvectors( this )
        !
        !< Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (SphericalUnitVectors), intent (in out)    :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_sphericalunitvectors
    !
    
    !
end module type_SphericalUnitVectors
