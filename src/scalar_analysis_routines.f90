module scalar_analysis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SpherepackUtility, only: &
        SpherepackUtility

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: shaec, shaeci
    public :: shaes, shaesi
    public :: shagc, shagci
    public :: shags, shagsi
    public :: ScalarAnalysisUtility

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    type, public :: ScalarAnalysisUtility
    contains
        ! Type-bound procedures
        procedure, nopass :: shaec
        procedure, nopass :: shaeci
        procedure, nopass :: shagc
        procedure, nopass :: shagci
        procedure, nopass :: shaes
        procedure, nopass :: shaesi
        procedure, nopass :: shags
        procedure, nopass :: shagsi
        procedure, nopass :: initialize_shaec
        procedure, nopass :: initialize_shaes
        procedure, nopass :: initialize_shagc
        procedure, nopass :: initialize_shags
    end type ScalarAnalysisUtility

    ! Declare interfaces for submodule implementation
    interface
        module subroutine shaec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshaec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshaec(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shaec

        module subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshaes, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshaes(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shaes

        module subroutine shagc(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshagc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshagc(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shagc

        module subroutine shags(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshags, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshags(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shags

        module subroutine shaeci(nlat, nlon, wshaec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshaec(:)
            integer(ip), intent(out) :: ierror
        end subroutine shaeci

        module subroutine shaesi(nlat, nlon, wshaes, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshaes(:)
            integer(ip), intent(out) :: ierror
        end subroutine shaesi

        module subroutine shagci(nlat, nlon, wshagc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshagc(:)
            integer(ip), intent(out) :: ierror
        end subroutine shagci

        module subroutine shagsi(nlat, nlon, wshags, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            real(wp),    intent(out)  :: wshags(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shagsi
    end interface

contains

    subroutine initialize_shaec(nlat, nlon, wshaec, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshaec(:)
        integer(ip),           intent(out) :: error_flag
        ! Local variables
        integer(ip) :: lshaec

        ! Get required workspace size
        lshaec = get_lshaec(nlat, nlon)

        ! Allocate memory
        allocate (wshaec(lshaec))

        ! Initialize wavetable
        call shaeci(nlat, nlon, wshaec, error_flag)

    end subroutine initialize_shaec

    subroutine initialize_shaes(nlat, nlon, wshaes, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshaes(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lshaes

        ! Get required workspace size
        lshaes = get_lshaes(nlat, nlon)

        ! Allocate memory
        allocate (wshaes(lshaes))

        ! Initialize wavetable
        call shaesi(nlat, nlon, wshaes, error_flag)

    end subroutine initialize_shaes

    subroutine initialize_shagc(nlat, nlon, wshagc, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshagc(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lshagc

        ! Get required workspace size
        lshagc = get_lshagc(nlat, nlon)

        ! Allocate memory
        allocate (wshagc(lshagc))

        ! Initialize wavetable
        call shagci(nlat, nlon, wshagc, error_flag)

    end subroutine initialize_shagc

    subroutine initialize_shags(nlat, nlon, wshags, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshags(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lshags

        ! Get required workspace size
        lshags = get_lshags(nlat, nlon)

        ! Allocate memory
        allocate (wshags(lshags))

        ! Initialize wavetable
        call shagsi(nlat, nlon, wshags, error_flag)

    end subroutine initialize_shags

    pure function get_lshaec(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip)             :: n1, n2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, n1, n2)

        return_value = 2*nlat*n2+3*((n1-2)*(2*nlat-n1-1))/2+nlon+15

    end function get_lshaec

    pure function get_lshagc(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip)             :: n1, n2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, n1, n2)

        return_value = nlat*(2*n2+3*n1-2)+3*n1*(1-n1)/2+nlon+15

    end function get_lshagc

    pure function get_lshaes(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip)             :: n1, n2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, n1, n2)

        return_value = (n1 * n2 * (2*nlat-n1+1))/2 + (nlon + 15)

    end function get_lshaes

    pure function get_lshags(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip)             :: n1, n2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, n1, n2)

        return_value = nlat*(3*(n1+n2)-2)+(n1-1)*(n2*(2*nlat-n1)-3*n1)/2+nlon+15

    end function get_lshags

end module scalar_analysis_routines
