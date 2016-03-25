module type_Grid

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: SphericalGrid

    !----------------------------------------------------------------------
    ! Dictionary: global variables confined to the module
    !----------------------------------------------------------------------
    integer (ip) :: allocate_status !! To check allocation status
    integer (ip) :: deallocate_status !! To check deallocation status
    !----------------------------------------------------------------------

    ! Declare derived data type
    type, abstract, public ::  SphericalGrid
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


    subroutine destroy_grid( this )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid), intent (in out) :: this
        !----------------------------------------------------------------------

        if ( this%initialized .eqv. .false. ) return

        ! Release memory
        if (allocated(this%grid_type)) deallocate(this%grid_type)
        if (allocated(this%longitudes)) deallocate(this%longitudes)
        if (allocated(this%latitudes)) deallocate(this%latitudes)

        ! Reset constants
        this%NUMBER_OF_LONGITUDES = 0
        this%NUMBER_OF_LATITUDES = 0
        this%LONGITUDINAL_MESH = 0.0_wp

        ! Reset flag
        this%initialized = .false.

    end subroutine destroy_grid


    subroutine get_equally_spaced_longitudes( this, nlon, phi )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid),           intent (in out) :: this
        integer (ip),           intent (in)     :: nlon !! number of longitudinal points
        real (wp), allocatable, intent (out)    :: phi(:)  !! longitudes: 0 <= phi <= 2*pi
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)         :: l !! counter
        real (wp), parameter :: TWO_PI = 2.0_wp * acos( -1.0_wp )
        !----------------------------------------------------------------------

        if ( nlon <= 0 ) then
            error stop 'TYPE (Grid): '&
                //'invalid argument NLON in GET_EQUALLY_SPACED_LONGITUDES'
        end if

        ! Release memory
        if (allocated(phi)) deallocate( phi, stat=deallocate_status )
        ! Check allocation status
        if ( deallocate_status /= 0 ) then
            error stop 'TYPE (Grid): '&
                //'Deallocating PHI failed in GET_EQUALLY_SPACED_LONGITUDES'
        end if

        ! Allocate memory
        allocate( phi(nlon), stat=allocate_status )
        ! Check allocation status
        if ( allocate_status /= 0 ) then
            error stop 'TYPE (Grid): '&
                //'Allocating PHI failed in GET_EQUALLY_SPACED_LONGITUDES'
        end if

        associate( dphi => this%LONGITUDINAL_MESH )
            ! Set equally spaced (uniform) mesh size
            dphi= TWO_PI / nlon
            ! Compute  equally spaced (uniform) longitudinal grid
            do l = 1, nlon
                phi(l) = real(l - 1, kind=wp) * dphi
            end do
        end associate

    end subroutine get_equally_spaced_longitudes


    subroutine print_to_unformatted_binary_files( this, header )
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (SphericalGrid),       intent (in out) :: this
        character (len=*),  intent (in)     :: header
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)  :: file_unit
        integer (ip)  :: record_length
        !----------------------------------------------------------------------

        ! Check if object is usable
        if ( this%initialized .eqv. .false. ) then
            error stop 'TYPE(Grid): '&
                //'uninitialized object in PRINT_TO_UNFORMATTED_BINARY_FILES'
        end if

        ! Write latitudes
        associate( theta => this%latitudes )
            inquire( iolength=record_length ) theta
            open( newunit=file_unit, file=header//this%grid_type//'_latitudes.dat', status='replace', &
                form='unformatted', access='direct', recl=record_length )
            write( file_unit, rec=1 ) theta
            close( file_unit )
        end associate

        ! Write longitudes
        associate( phi => this%longitudes )
            inquire( iolength=record_length ) phi
            open( newunit=file_unit, file=header//this%grid_type//'_longitudes.dat', status='replace', &
                form='unformatted', access='direct', recl=record_length )
            write( file_unit, rec=1 ) phi
            close( file_unit )
        end associate

    end subroutine print_to_unformatted_binary_files

end module type_Grid
