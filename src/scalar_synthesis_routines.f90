module scalar_synthesis_routines

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
    public :: shses, shsesi
    public :: shsec, shseci
    public :: shsgs, shsgsi
    public :: shsgc, shsgci
    public :: ScalarSynthesisUtility

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    type, public :: ScalarSynthesisUtility
    contains
        ! Type-bound procedures
        procedure, nopass :: shsec
        procedure, nopass :: shseci
        procedure, nopass :: shsgc
        procedure, nopass :: shsgci
        procedure, nopass :: shses
        procedure, nopass :: shsesi
        procedure, nopass :: shsgs
        procedure, nopass :: shsgsi
        procedure, nopass :: initialize_shsec
        procedure, nopass :: initialize_shses
        procedure, nopass :: initialize_shsgc
        procedure, nopass :: initialize_shsgs
    end type ScalarSynthesisUtility

    ! Declare interfaces for submodule implementation
    interface
        module subroutine shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshses(:)
            integer(ip), intent(out) :: ierror
        end subroutine shses

        module subroutine shsesi(nlat, nlon, wshses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshses(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsesi

        module subroutine shsgs(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, wshsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: mode
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsgs

        module subroutine shsgsi(nlat, nlon, wshsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            real(wp),    intent(out)  :: wshsgs(:)
            integer(ip), intent(out)  :: ierror
        end subroutine shsgsi

        module subroutine shsec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsec

        module subroutine shseci(nlat, nlon, wshsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine shseci

        module subroutine shsgc(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, wshsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: mode
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsgc

        module subroutine shsgci(nlat, nlon, wshsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine shsgci
    end interface

contains

    subroutine initialize_shsec(nlat, nlon, wshsec, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshsec(:)
        integer(ip),           intent(out) :: error_flag
        ! Local variables
        integer(ip) :: lshsec

        ! Get required workspace size
        lshsec = get_lshsec(nlat, nlon)

        ! Allocate memory
        allocate (wshsec(lshsec))

        ! Initialize wavetable
        call shseci(nlat, nlon, wshsec, error_flag)

    end subroutine initialize_shsec

    subroutine initialize_shses(nlat, nlon, wshses, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshses(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lshses

        ! Get required workspace size
        lshses = get_lshses(nlat, nlon)

        ! Allocate memory
        allocate (wshses(lshses))

        ! Initialize wavetable
        call shsesi(nlat, nlon, wshses, error_flag)

    end subroutine initialize_shses

    subroutine initialize_shsgc(nlat, nlon, wshsgc, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshsgc(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lshsgc

        ! Get required workspace size
        lshsgc = get_lshsgc(nlat, nlon)

        ! Allocate memory
        allocate (wshsgc(lshsgc))

        ! Initialize wavetable
        call shsgci(nlat, nlon, wshsgc, error_flag)

    end subroutine initialize_shsgc

    subroutine initialize_shsgs(nlat, nlon, wshsgs, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wshsgs(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lshsgs

        ! Get required workspace size
        lshsgs = get_lshsgs(nlat, nlon)

        ! Allocate memory
        allocate (wshsgs(lshsgs))

        ! Initialize wavetable
        call shsgsi(nlat, nlon, wshsgs, error_flag)

    end subroutine initialize_shsgs

    pure function get_lshsec(nlat, nlon) &
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

    end function get_lshsec

    pure function get_lshsgc(nlat, nlon) &
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

    end function get_lshsgc

    pure function get_lshses(nlat, nlon) &
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

    end function get_lshses

    pure function get_lshsgs(nlat, nlon) &
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

    end function get_lshsgs

end module scalar_synthesis_routines
