module vector_synthesis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_SpherepackUtility, only: &
        SpherepackUtility

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vhsgc, vhsgci, initialize_vhsec
    public :: vhses, vhsesi, initialize_vhses
    public :: vhsec, vhseci, initialize_vhsgc
    public :: vhsgs, vhsgsi, initialize_vhsgs
    public :: VectorSynthesisUtility
    public :: get_lvhsgs
    
    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    type, public :: VectorSynthesisUtility
    contains
        ! Type-bound procedures
        procedure, nopass :: vhsec
        procedure, nopass :: vhseci
        procedure, nopass :: vhsgc
        procedure, nopass :: vhsgci
        procedure, nopass :: vhses
        procedure, nopass :: vhsesi
        procedure, nopass :: vhsgs
        procedure, nopass :: vhsgsi
        procedure, nopass :: initialize_vhsec
        procedure, nopass :: initialize_vhses
        procedure, nopass :: initialize_vhsgc
        procedure, nopass :: initialize_vhsgs
    end type VectorSynthesisUtility

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vhsec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsec

        module subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhses(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhses

        module subroutine vhsgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgc

        module subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgs

        module subroutine vhseci(nlat, nlon, wvhsec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhsec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhseci

        module subroutine vhsesi(nlat, nlon, wvhses, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhses(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsesi

        module subroutine vhsgci(nlat, nlon, wvhsgc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhsgc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgci

        module subroutine vhsgsi(nlat, nlon, wvhsgs, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhsgs(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhsgsi
    end interface

contains

    subroutine initialize_vhsec(nlat, nlon, wvhsec, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wvhsec(:)
        integer(ip),           intent(out) :: error_flag
        ! Local variables
        integer(ip) :: lvhsec

        ! Get required workspace size
        lvhsec = get_lvhsec(nlat, nlon)

        ! Allocate memory
        allocate (wvhsec(lvhsec))

        ! Initialize wavetable
        call vhseci(nlat, nlon, wvhsec, error_flag)

    end subroutine initialize_vhsec

    subroutine initialize_vhses(nlat, nlon, wvhses, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wvhses(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lvhses

        ! Get required workspace size
        lvhses = get_lvhses(nlat, nlon)

        ! Allocate memory
        allocate (wvhses(lvhses))

        ! Initialize wavetable
        call vhsesi(nlat, nlon, wvhses, error_flag)

    end subroutine initialize_vhses

    subroutine initialize_vhsgc(nlat, nlon, wvhsgc, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wvhsgc(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lvhsgc

        ! Get required workspace size
        lvhsgc = get_lvhsgc(nlat, nlon)

        ! Allocate memory
        allocate (wvhsgc(lvhsgc))

        ! Initialize wavetable
        call vhsgci(nlat, nlon, wvhsgc, error_flag)

    end subroutine initialize_vhsgc

    subroutine initialize_vhsgs(nlat, nlon, wvhsgs, error_flag)

        ! Dummy arguments
        integer(ip),           intent(in)  :: nlat
        integer(ip),           intent(in)  :: nlon
        real(wp), allocatable, intent(out) :: wvhsgs(:)
        integer(ip),           intent(out) :: error_flag

        ! Local variables
        integer(ip) :: lvhsgs

        ! Get required workspace size
        lvhsgs = get_lvhsgs(nlat, nlon)

        ! Allocate memory
        allocate (wvhsgs(lvhsgs))

        ! Initialize wavetable
        call vhsgsi(nlat, nlon, wvhsgs, error_flag)

    end subroutine initialize_vhsgs

    pure function get_lvhsgc(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip)             :: n1, n2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, n1, n2)

        return_value = 4 * nlat * n2 + 3 * max(n1-2,0)*(2*nlat-n1-1) + nlon + 15

    end function get_lvhsgc

    pure function get_lvhsgs(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip)  :: imid, lmn

        imid = (nlat+1)/2
        lmn = (nlat*(nlat+1))/2
        return_value = 2*(imid*lmn)+nlon+15

    end function get_lvhsgs

    pure function get_lvhsec(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip)             :: n1, n2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, n1, n2)

        return_value = 4*nlat*n2+3*max(n1-2, 0)*(nlat+nlat-n1-1)+nlon+15

    end function get_lvhsec

    pure function get_lvhses(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip)             :: n1, n2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, n1, n2)

        return_value = n1 * n2 * ((2*nlat) - n1 + 1) + nlon + 15

    end function get_lvhses

end module vector_synthesis_routines
