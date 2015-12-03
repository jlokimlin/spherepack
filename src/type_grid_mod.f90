module type_grid_mod

    use, intrinsic :: iso_fortran_env, only: &
        REAL64, &
        INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private

    ! Public derived data type
    public :: grid_t

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    integer, parameter    :: WP = REAL64          !! 64 bit real
    integer, parameter    :: IP = INT32           !! 32 bit integer
    character(200)        :: error_message        !! Probably long enough
    integer (IP)          :: allocate_status      !! To check allocation status
    integer (IP)          :: deallocate_status    !! To check deallocation status
    real (WP), parameter  :: TWO_PI = 8.0_WP * atan( 1.0_WP )
    real (WP), parameter  :: PI     = 4.0_WP * atan( 1.0_WP )
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type ::  grid_t

        ! Components
        logical, private                      :: initialized      = .false. !! Flag to check if object is instantiated

        character(3)                          :: grid_type        = 'GAU'   !! either 'REG' or 'GAU'

        real (WP)                             :: mesh_phi         = 0.0_WP  !! Uniform mesh in phi
        real (WP)                             :: mesh_theta       = 0.0_WP  !! Only used in 'REG' grid

        real (WP), dimension (:), allocatable :: gaussian_weights           !! Used for integration, requires 'GAU' for allocation
        real (WP), dimension (:), allocatable :: latitudes                  !! 0 <= theta <= pi
        real (WP), dimension (:), allocatable :: longitudes                 !! 0 <= phi <= 2*pi


    contains

        ! All methods are private unless stated otherwise
        private

        ! Public methods
        procedure, non_overridable, public :: Create
        procedure, non_overridable, public :: Destroy
        procedure, nopass, public          :: Get_gaussian_weights_and_points !! Gaqd

        ! Private methods
        procedure                          :: Get_equally_spaced_longitudes
        procedure                          :: Get_equally_spaced_latitudes
        final                              :: Finalize

    end type grid_t

contains
    !
    !*****************************************************************************************
    !
    subroutine Create( this, nlat, nlon, grid_type )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), intent (in out)     :: this
        integer (IP), intent (in)           :: nlat         !! number of latitudinal points
        integer (IP), intent (in)           :: nlon         !! number of longitudinal points
        character(*), intent (in), optional :: grid_type    !! Either 'GAU' or 'REG'
        !--------------------------------------------------------------------------------

        ! Check status
        if ( this%initialized ) then
            print *, 'ERROR: You must destroy "grid" before re-instantiating'
            return
        end if

        ! Set latitudes: 0 <= theta <= pi
        if ( present(grid_type) ) then

            ! Check if the grid type is regular
            if ( grid_type .eq. 'REG') then

                ! Set the grid type
                this%grid_type = grid_type

                ! Compute latitudes: 0 <= theta <= pi
                call this%Get_equally_spaced_latitudes( &
                    nlat, this%latitudes )

            ! Check if the grid type is gaussian
            else if ( grid_type .eq. 'GAU') then

                call this%Get_gaussian_weights_and_points( &
                    nlat, this%latitudes, this%gaussian_weights )

            else

                ! Handle invalid grid type
                print *, 'ERROR: optional argument grid_type = ', grid_type
                print *, 'must be either GAU or REG (default GAU)'
                stop

            end if

        else

            ! Set default '(GAU)' grid
            call this%Get_gaussian_weights_and_points( &
                nlat, this%latitudes, this%gaussian_weights )

        end if

        ! Compute longitudes: 0 <= phi <= 2*pi
        call this%Get_equally_spaced_longitudes( nlon, this%longitudes )

        ! Set status
        this%initialized = .true.

    end subroutine Create
    !
    !*****************************************************************************************
    !
    subroutine Destroy( this )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        ! Check status
        if ( .not. this%initialized ) return

        ! Reset default grid type
        this%grid_type = 'GAU'

        ! Reset meshes
        this%mesh_phi   = 0.0_WP
        this%mesh_theta = 0.0_WP

        ! Destroy gaussian_weights
        if ( allocated( this%gaussian_weights ) ) then

            deallocate ( &
                this%gaussian_weights, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed '&
                    &//'for "gaussian_weights":', &
                    trim( error_message )
                return
            end if
        end if

        ! Destroy latitudes
        if ( allocated( this%latitudes ) ) then

            deallocate ( &
                this%latitudes, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed '&
                    &//'for "latitudes":', &
                    trim( error_message )
                return
            end if
        end if

        ! Destroy longitudes
        if ( allocated( this%longitudes ) ) then
            deallocate ( &
                this%longitudes, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocation failed '&
                    &//'for "longitudes":', &
                    trim( error_message )
                return
            end if
        end if

        ! Reset status
        this%initialized = .false.

    end subroutine Destroy
    !
    !*****************************************************************************************
    !
    subroutine Get_gaussian_weights_and_points( nlat, theta, wts )
        !
        ! Purpose:
        !
        ! Computes the nlat-many gaussian (co)latitudes and weights.
        ! the colatitudes are in radians and lie in the interval (0,pi).
        !
        ! References:
        !
        ! [1] Swarztrauber, Paul N.
        !     "On computing the points and weights for Gauss--Legendre quadrature."
        !     SIAM Journal on Scientific Computing 24.3 (2003): 945-954.

        ! [2]  http://www2.cisl.ucar.edu/spherepack/documentation#gaqd.html
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        integer (IP), intent (in)             :: nlat  !! number of latitudinal points
        real (WP), dimension (:), allocatable :: theta !! latitudinal points: 0 <= theta <= pi
        real (WP), dimension (:), allocatable :: wts   !! gaussian weights
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)  :: ierror
        integer (IP)  :: dummy_integer !! unused integer variable to maintain backwards compatibility
        real (WP)     :: dummy_real    !! unused double precision variable to maintain backwards compatibility
        !--------------------------------------------------------------------------------

        ! Allocate latitudes
        if ( allocated( theta ) ) then

            ! Deallocate if necessary
            deallocate ( &
                theta, &
                stat = deallocate_status, &
                errmsg = error_message )

            ! Check deallocate status
            if ( deallocate_status /= 0 ) then
                print *, 'Deallocating theta failed in '&
                    &//'in "Get_gaussian_weights_and_points"', &
                    trim( error_message )
                return
            end if

        else

        allocate ( &
            theta( 1:nlat ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            print *, 'Allocating theta failed in '&
                &//'in "Get_gaussian_weights_and_points"', &
                trim( error_message )
            return
        end if

    end if

    ! Allocate gaussian weights
    if ( allocated( wts ) ) then

        ! Deallocate if necessary
        deallocate ( &
            wts, &
            stat = deallocate_status, &
            errmsg = error_message )

        ! Check deallocate status
        if ( deallocate_status /= 0 ) then
            print *, 'Deallocating wts failed in '&
                &//'in "Get_gaussian_weights_and_points"', &
                trim( error_message )
            return
        end if

    else

    allocate ( &
        wts( 1:nlat ), &
        stat = allocate_status, &
        errmsg = error_message )

    ! Check allocate status
    if ( allocate_status /= 0 ) then
        print *, 'Allocating wts failed in '&
            &//'in "Get_gaussian_weights_and_points"', &
            trim( error_message )
        return
    end if

end if

!dummy_integer = nlat * (nlat + 2)

! Set latitudes and gaussian weights
call Gaqd( nlat, theta, wts, dummy_real, dummy_integer, ierror )

! Check error flag
if ( ierror /= 0 ) then
    print *, 'SPHEREPACK 3.2 error = ', ierror,' in Gaqd'
end if

end subroutine Get_gaussian_weights_and_points
!
!*****************************************************************************************
!
subroutine Get_equally_spaced_latitudes( this, nlat, theta )
    !
    ! Purpose:
    !--------------------------------------------------------------------------------
    ! Dictionary: calling arguments
    !--------------------------------------------------------------------------------
    class (grid_t), intent (in out)         :: this
    integer (IP), intent (in)               :: nlat  !! number of latitudinal points
    real (WP), dimension (:), allocatable   :: theta !! latitudes: 0 <= theta <= pi
    !--------------------------------------------------------------------------------
    ! Dictionary: local variables
    !--------------------------------------------------------------------------------
    integer (IP) :: k     !! counter
    !--------------------------------------------------------------------------------

    ! Allocate array
    if ( allocated( theta ) ) then

        ! Deallocate if already allocated
        deallocate ( &
            theta, &
            stat = deallocate_status, &
            errmsg = error_message )

        ! Check deallocation status
        if ( deallocate_status /= 0 ) then
            print *, 'Deallocation failed '&
                &//'in Get_equally_spaced_latitudes: ', &
                trim( error_message )
            return
        end if

    else

    ! Allocate array
    allocate ( &
        theta( 1:nlat ), &
        stat = allocate_status, &
        errmsg = error_message )

    ! Check allocation status
    if ( allocate_status /= 0 ) then
        print *, 'Allocation failed '&
            &//'in Get_equally_spaced_latitudes: ', &
            trim( error_message )
        return
    end if

end if

! Set equally spaced (uniform) mesh size
this%mesh_theta = PI / nlat

associate( Dtheta => this%mesh_theta )

    ! Compute latitudinal grid
    do concurrent ( k = 1:nlat )

        theta(k) = real(k - 1, WP) * Dtheta

    end do

end associate

end subroutine Get_equally_spaced_latitudes
!
!*****************************************************************************************
!
subroutine Get_equally_spaced_longitudes( this, nlon, phi )
    !
    ! Purpose:
    !--------------------------------------------------------------------------------
    ! Dictionary: calling arguments
    !--------------------------------------------------------------------------------
    class (grid_t), intent (in out)         :: this
    integer (IP), intent (in)               :: nlon !! number of longitudinal points
    real (WP), dimension (:), allocatable   :: phi  !! longitudes: 0 <= phi <= 2*pi
    !--------------------------------------------------------------------------------
    ! Dictionary: local variables
    !--------------------------------------------------------------------------------
    integer (IP)  :: l   !! counter
    !--------------------------------------------------------------------------------

    ! Allocate array
    if ( allocated( phi ) ) then

        ! Deallocate if already allocated
        deallocate ( &
            phi, &
            stat = deallocate_status, &
            errmsg = error_message )

        ! Check deallocation status
        if ( deallocate_status /= 0 ) then
            print *, 'Deallocation failed '&
                &//'in Get_equally_spaced_longitudes: ', &
                trim( error_message )
            return
        end if

    else

    ! Allocate array
    allocate ( &
        phi( 1:nlon ), &
        stat = allocate_status, &
        errmsg = error_message )

    ! Check allocation status
    if ( allocate_status /= 0 ) then
        print *, 'Allocation failed '&
            &//'in Get_equally_spaced_longitudes: ', &
            trim( error_message )
        return
    end if

end if

! Set equally spaced (uniform) mesh size
this%mesh_phi = TWO_PI / nlon

associate( Dphi => this%mesh_phi )

    ! Compute longitudinal grid
    do concurrent ( l = 1:nlon )

        phi(l) = real(l - 1, WP) * Dphi

    end do

end associate

end subroutine Get_equally_spaced_longitudes
!
!*****************************************************************************************
!
subroutine Finalize( this )
    !
    ! Purpose:
    !< Finalize object
    !
    !--------------------------------------------------------------------------------
    ! Dictionary: calling arguments
    !--------------------------------------------------------------------------------
    type (grid_t), intent (in out) :: this
    !--------------------------------------------------------------------------------

    call this%Destroy()

end subroutine Finalize
    !
    !*****************************************************************************************
    !
end module type_grid_mod
