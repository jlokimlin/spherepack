module type_grid_mod

    use, intrinsic :: iso_fortran_env, only: &
        REAL64, &
        INT32

    ! Explicit typing only
    implicit none

    ! Everything is private except the derived data type itself
    private

    !---------------------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !---------------------------------------------------------------------------------
    integer, parameter    :: WP     = REAL64                    !! 64 bit real
    integer, parameter    :: IP     = INT32                     !! 32 bit integer
    real (WP)             :: PI     = 4.0_WP * atan( 1.0_WP )   !! To calculate latitudes:   0 <= theta <= pi
    real (WP)             :: TWO_PI = 8.0_WP * atan( 1.0_WP )   !! To calculate longitudes:  0 <=  phi  <= 2*pi
    character(200)        :: error_message                      !! Probably long enough
    integer (IP)          :: allocate_status                    !! To check allocation status
    integer (IP)          :: deallocate_status                  !! To check deallocation status
    !---------------------------------------------------------------------------------

    ! Declare derived data type
    type, public ::  grid_t

        ! Derived data type components
        logical, private                      :: initialized      = .false. !! Flag to check if object is instantiated
        character(3)                          :: grid_type        = 'GAU'   !! either 'REG' or 'GAU'
        real (WP)                             :: mesh_phi         = 0.0_WP  !! Uniform mesh in phi
        real (WP)                             :: mesh_theta       = 0.0_WP  !! Only used in 'REG' grid
        real (WP), dimension (:), allocatable :: gaussian_weights           !! Used for integration, requires 'GAU' for allocation
        real (WP), dimension (:), allocatable :: latitudes                  !! 0 <= theta <= pi
        real (WP), dimension (:), allocatable :: longitudes                 !! 0 <= phi <= 2*pi


    contains

        ! Public methods
        procedure, non_overridable :: Create
        procedure, non_overridable :: Destroy
        procedure                  :: Set_equally_spaced_longitudes
        procedure                  :: Set_equally_spaced_latitudes
        procedure                  :: Gaussian_weights_and_points !! Known as Gaqd in SPHEREPACK 3.2
        !final, non_overridable               :: Finalize

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
        integer (IP), intent (in)           :: nlat, nlon
        character(*), intent (in), optional :: grid_type
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

                this%grid_type = grid_type
                call this%Set_equally_spaced_latitudes( nlat )

            ! Check if the grid type is gaussian
            else if ( grid_type .eq. 'GAU') then

                call this%Gaussian_weights_and_points( nlat )

            else

                ! Handle invalid grid type
                print *, 'ERROR: optional argument grid_type = ', grid_type
                print *, 'must be either GAU or REG (default GAU)'
                stop

            end if

        else

            ! Set default '(GAU)' grid
            call this%Gaussian_weights_and_points( nlat )

        end if


        ! Set longitudes: 0 <= phi <= 2*pi
        call this%Set_equally_spaced_longitudes( nlon )

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
    subroutine Gaussian_weights_and_points( this, nlat )
        !
        ! References:
        !
        ! [1] Swarztrauber, Paul N.
        !     "On computing the points and weights for Gauss--Legendre quadrature."
        !     SIAM Journal on Scientific Computing 24.3 (2003): 945-954.

        ! [2]  https://www2.cisl.ucar.edu/gaussian_gridpack/documentation#gaqd.html
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), intent (in out)    :: this
        integer (IP), intent (in)          :: nlat
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)    :: LWORK, ierror
        real (WP)       :: dummy_variable ! dummy_variable is an unused
        ! double precision variable that permits a simple exchange with
        ! the old FORTRAN 77 routine with the same name in SPHEREPACK 3.2.
        !--------------------------------------------------------------------------------

        ! Check status
        if ( this%initialized ) then
            print *, 'ERROR: You must destroy "grid" before re-instantiating'
            return
        end if

        ! Allocate arrays
        allocate ( &
            this%latitudes( 1:nlat ), &
            this%gaussian_weights( 1:nlat ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed '&
                &//'in Gaussian_weights_and_points: ', &
                trim( error_message )
            return
        end if

        ! Set LWORK
        LWORK = nlat * (nlat + 2)

        ! Set latitudes and gaussian weights
        call Gaqd( &
            nlat, this%latitudes, &
            this%gaussian_weights, dummy_variable, &
            LWORK, ierror )

        ! Check error flag
        if ( ierror /= 0 ) then
            print *, 'SPHEREPACK 3.2 error = ', ierror,' in Gaqd'
            return
        end if

    end subroutine Gaussian_weights_and_points
    !
    !*****************************************************************************************
    !
    subroutine Set_equally_spaced_latitudes( this, nlat )
        !
        ! Purpose:
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), target, intent (in out) :: this
        integer (IP), intent (in)               :: nlat
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                      :: k               !! counter
        real (WP)                         :: mesh_theta      !! Uniform mesh in theta
        real (WP), dimension (:), pointer :: theta => null() !! address of latitudes
        !--------------------------------------------------------------------------------

        ! Check status
        if ( this%initialized ) then
            print *, 'ERROR: You must destroy "grid" before re-instantiating'
            return
        end if

        ! Set the uniform latitudinal mesh
        mesh_theta = PI / nlat
        this%mesh_theta = mesh_theta

        ! Allocate array
        allocate ( &
            this%latitudes( 1:nlat ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocate status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed '&
                &//'in Set_equally_spaced_latitudes: ', &
                trim( error_message )
            return
        end if

        ! Associate pointer
        if ( .not. associated( theta ) ) theta => this%latitudes

        ! Compute longitudinal grid
        do k = 1, size( theta )

            theta(k) = real(k - 1, WP) * mesh_theta

        end do

        ! Nullify pointer
        if ( associated(theta) ) nullify (theta)

    end subroutine Set_equally_spaced_latitudes
    !
    !*****************************************************************************************
    !
    subroutine Set_equally_spaced_longitudes( this, nlon )
        !
        ! Purpose:
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (grid_t), target, intent (in out) :: this
        integer (IP), intent (in)               :: nlon !! number of longitudes
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (IP)                      :: l             !! counter
        real (WP)                         :: mesh_phi      !! equally space (uniform) mesh
        real (WP), dimension (:), pointer :: phi => null() !! address of longitudes
        !--------------------------------------------------------------------------------

        ! Check status
        if ( this%initialized ) then
            print *, 'ERROR: You must destroy "grid" before re-instantiating'
            return
        end if

        ! Set equally spaced (uniform) mesh size
        mesh_phi = TWO_PI / real(nlon, WP)
        this%mesh_phi = mesh_phi

        ! Allocate array
        allocate ( &
            this%longitudes( 1:nlon ), &
            stat = allocate_status, &
            errmsg = error_message )

        ! Check allocation status
        if ( allocate_status /= 0 ) then
            print *, 'Allocation failed '&
                &//'in Set_longitudes: ', &
                trim( error_message )
            return
        end if

        ! Associate pointer
        if ( .not. associated( phi ) ) phi => this%longitudes

        ! Compute longitudinal grid
        do l = 1, size( phi )

            phi(l) = real(l - 1, WP) * mesh_phi

        end do

        ! Nullify pointer
        if ( associated( phi ) ) nullify ( phi )


    end subroutine Set_equally_spaced_longitudes
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
        class (grid_t), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%Destroy()

    end subroutine Finalize
    !
    !*****************************************************************************************
    !
end module type_grid_mod
