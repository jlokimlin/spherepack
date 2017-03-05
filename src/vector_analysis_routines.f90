module vector_analysis_routines

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
    public :: vhaec, vhaes, vhagc, vhags
    public :: vhaeci, vhaesi, vhagci, vhagsi

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhaec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhaec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaec

        module subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhaes, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhaes(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaes

        module subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhagc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(mdab, ndab, nt)
            real(wp),    intent(out) :: cr(mdab, ndab, nt)
            real(wp),    intent(out) :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhagc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhagc

        module subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhags, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: v(idvw, jdvw, nt)
            real(wp),    intent(in)  :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(out) :: br(mdab, ndab, nt)
            real(wp),    intent(out) :: bi(..)
            real(wp),    intent(out) :: cr(..)
            real(wp),    intent(out) :: ci(..)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvhags(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhags

        module subroutine vhaeci(nlat, nlon, wvhaec, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhaec(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaeci

        module subroutine vhaesi(nlat, nlon, wvhaes, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhaes(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhaesi

        module subroutine vhagci(nlat, nlon, wvhagc, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhagc(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhagci

        module subroutine vhagsi(nlat, nlon, wvhags, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhags(:)
            integer(ip), intent(out) :: ierror
        end subroutine vhagsi
    end interface

    type, public :: VhaesAux
    contains
        ! Type-bound procedures
        procedure, nopass :: vhaes
        procedure, nopass :: vhaesi
        procedure, nopass :: get_lvhaes
    end type VhaesAux

    type, public :: VhagsAux
    contains
        ! Type-bound procedures
        procedure, nopass :: vhags
        procedure, nopass :: vhagsi
        procedure, nopass :: get_lvhags
    end type VhagsAux

contains

    pure function get_lvhags(nlat, nlon) result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables
        integer(ip)         :: l1, l2
        type(SpherepackUtility) :: util

        call util%compute_parity(nlat, nlon, l1, l2)

        return_value = ((nlat+1)**2)*nlat/2+nlon+15
        !return_value = max(3*nlat*(nlat+1)+2, l1*l2*(2*nlat-l1+1)+nlon+15+2*nlat)

    end function get_lvhags

    pure function get_lvhaes(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackUtility) :: util


        call util%compute_parity(nlat, nlon, l1, l2)

        return_value = l1*l2*(2*nlat-l1+1)+nlon+15

    end function get_lvhaes

end module vector_analysis_routines
