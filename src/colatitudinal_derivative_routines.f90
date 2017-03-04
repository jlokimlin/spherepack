module colatitudinal_derivative_routines

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
    public :: vtsec, vtses, vtsgc, vtsgs
    public :: vtseci, vtsesi, vtsgci, vtsgsi

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vtsec(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, lwvts, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vtsec

        module subroutine vtses(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, lwvts, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vtses

        module subroutine vtsgc(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, lwvts, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vtsgc

        module subroutine vtsgs(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvts, lwvts, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vt(idvw, jdvw, nt)
            real(wp),    intent(out) :: wt(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdab, ndab, nt)
            real(wp),    intent(in)  :: bi(mdab, ndab, nt)
            real(wp),    intent(in)  :: cr(mdab, ndab, nt)
            real(wp),    intent(in)  :: ci(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vtsgs

        module subroutine vtseci(nlat, nlon, wvts, lwvts, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vtseci

        module subroutine vtsesi(nlat, nlon, wvts, lwvts, work, lwork, dwork, ldwork, &
            ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vtsesi

        module subroutine vtsgci(nlat, nlon, wvts, lwvts, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vtsgci

        module subroutine vtsgsi(nlat, nlon, wvts, lwvts, work, lwork, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvts(lwvts)
            integer(ip), intent(in)  :: lwvts
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vtsgsi
    end interface

end module colatitudinal_derivative_routines
