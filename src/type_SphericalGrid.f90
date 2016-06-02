module type_SphericalGrid

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: SphericalGrid


    ! Declare derived data type
    type, abstract, public :: SphericalGrid
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        logical,                        public :: initialized = .false.
        integer (ip),                   public :: NUMBER_OF_LONGITUDES = 0  !! number of longitudinal points
        integer (ip),                   public :: NUMBER_OF_LATITUDES = 0 !! number of latitudinal points
        real (wp),                      public :: LONGITUDINAL_MESH = 0.0_wp !! Only used in 'REG' grid
        real (wp),         allocatable, public :: latitudes(:)  !! 0 <= theta <= pi
        real (wp),         allocatable, public :: longitudes(:) !! 0 <= phi <= 2*p
        character (len=:), allocatable, public :: grid_type
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, public :: destroy_grid
        procedure, public :: get_equally_spaced_longitudes
        procedure, public :: print_to_unformatted_binary_files
        !----------------------------------------------------------------------
    end type SphericalGrid



contains



    subroutine destroy_grid(this)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        ! Check flag
        if (this%initialized .eqv. .false.) then
            return
        end if

        !
        !==> Release memory
        !
        if (allocated(this%grid_type)) then
            deallocate( this%grid_type )
        end if

        if (allocated(this%longitudes)) then
            deallocate( this%longitudes )
        end if

        if (allocated(this%latitudes)) then
            deallocate( this%latitudes )
        end if

        ! Reset constants
        this%NUMBER_OF_LONGITUDES = 0
        this%NUMBER_OF_LATITUDES = 0
        this%LONGITUDINAL_MESH = 0.0_wp

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_grid



    subroutine get_equally_spaced_longitudes(this, nlon, phi)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid),  intent (in out) :: this
        integer (ip),           intent (in)     :: nlon !! number of longitudinal points
        real (wp), allocatable, intent (out)    :: phi(:)  !! longitudes: 0 <= phi <= 2*pi
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)         :: i !! counter
        real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
        !----------------------------------------------------------------------

        !
        !==> Check validity of calling argument
        !
        if (nlon <= 0) then
            error stop 'Object of class (SphericalGrid): '&
                //'invalid calling argument nlon <= 0 '&
                //'in get_equally_spaced_longitudes'
        end if

        !
        !==> Allocate memory
        !
        allocate( phi(nlon) )

        !
        !==> Compute equally space (uniform) longitudinal grid
        !
        associate( dphi => this%LONGITUDINAL_MESH )

            ! Set equally spaced (uniform) mesh size
            dphi= TWO_PI / nlon

            ! Compute grid
            phi = [ (real(i - 1, kind=wp) * dphi, i=1, nlon) ]

        end associate

    end subroutine get_equally_spaced_longitudes



    subroutine print_to_unformatted_binary_files(this, header)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid), intent (in out) :: this
        character (len=*),     intent (in)     :: header
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: file_unit
        !----------------------------------------------------------------------

        ! Check if object is usable
        if (this%initialized .eqv. .false.) then
            error stop 'Uninitialized object of class (SphericalGrid): '&
                //'in print_to_unformatted_binary_files'
        end if

        ! Write latitudes
        associate( theta => this%latitudes )

            open( newunit=file_unit, &
                file=header//this%grid_type//'_latitudes.dat', &
                status='replace', action='write', &
                form='unformatted', access='stream' )
            write( file_unit ) theta
            close( file_unit )

        end associate

        ! Write longitudes
        associate( phi => this%longitudes )

            open( newunit=file_unit, &
                file=header//this%grid_type//'_longitudes.dat', &
                status='replace', action='write', &
                form='unformatted', access='stream' )
            write( file_unit ) phi
            close( file_unit )

        end associate

    end subroutine print_to_unformatted_binary_files


end module type_SphericalGrid
