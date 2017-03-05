module vector_laplacian_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use vector_synthesis_routines, only: &
        vhses, vhsec, vhsgc, vhsgs

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vlapec, vlapes, vlapgc, vlapgs
    public :: ivlapec, ivlapes, ivlapgc, ivlapgs
    public :: get_workspace_indices_for_inversion, perform_setup_for_inversion
    public :: perform_setup_for_vector_laplacian

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vlapec(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhsec, lvhsec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsec(lvhsec)
            integer(ip), intent(in)  :: lvhsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vlapec

        module subroutine vlapes(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhses, lvhses, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhses(lvhses)
            integer(ip), intent(in)  :: lvhses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vlapes

        module subroutine vlapgc(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhsgc, lvhsgc, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgc(lvhsgc)
            integer(ip), intent(in)  :: lvhsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vlapgc

        module subroutine vlapgs(nlat, nlon, ityp, nt, vlap, wlap, idvw, jdvw, br, bi, &
            cr, ci, mdbc, ndbc, wvhsgs, lvhsgs, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vlap(idvw, jdvw, nt)
            real(wp),    intent(out) :: wlap(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgs(lvhsgs)
            integer(ip), intent(in)  :: lvhsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vlapgs

        module subroutine ivlapec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhsec, lvhsec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsec(lvhsec)
            integer(ip), intent(in)  :: lvhsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine ivlapec

        module subroutine ivlapes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhses, lvhses, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhses(lvhses)
            integer(ip), intent(in)  :: lvhses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine ivlapes

        module subroutine ivlapgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhsgc, lvhsgc, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgc(lvhsgc)
            integer(ip), intent(in)  :: lvhsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine ivlapgc

        module subroutine ivlapgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdbc, ndbc, wvhsgs, lvhsgs, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: ityp
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: br(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: bi(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: cr(mdbc, ndbc, nt)
            real(wp),    intent(in)  :: ci(mdbc, ndbc, nt)
            integer(ip), intent(in)  :: mdbc
            integer(ip), intent(in)  :: ndbc
            real(wp),    intent(in)  :: wvhsgs(lvhsgs)
            integer(ip), intent(in)  :: lvhsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine ivlapgs
    end interface

contains

    pure subroutine compute_coefficient_multipliers(fnn)

        ! Dummy arguments
        real(wp), intent(out) :: fnn(:)

        ! Local variables
        integer(ip) :: n
        real(wp)    :: fn

        associate (nlat => size(fnn))
            do n=2, nlat
                fn = real(n - 1, kind=wp)
                fnn(n) = -fn * (fn + ONE)
            end do
        end associate

    end subroutine compute_coefficient_multipliers

    pure subroutine perform_setup_for_vector_laplacian(ityp, brlap, bilap, crlap, cilap, br, bi, cr, ci, fnn)

        ! Dummy arguments
        integer(ip), intent(out) :: ityp
        real(wp),    intent(in)  :: br(:, :, :)
        real(wp),    intent(in)  :: bi(:, :, :)
        real(wp),    intent(in)  :: cr(:, :, :)
        real(wp),    intent(in)  :: ci(:, :, :)
        real(wp),    intent(out) :: brlap(:, :, :)
        real(wp),    intent(out) :: bilap(:, :, :)
        real(wp),    intent(out) :: crlap(:, :, :)
        real(wp),    intent(out) :: cilap(:, :, :)
        real(wp),    intent(out) :: fnn(:)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            mmax => size(brlap, dim=1), &
            nlat => size(brlap, dim=2), &
            nt => size(brlap, dim=3) &
           )

            ! Preset coefficient multiplyers in vector
            call compute_coefficient_multipliers(fnn)

            !  Set vector laplacian coefficients from br, bi, cr, ci
            select case (ityp)
                case (0, 3, 6)
                    !
                    !     all coefficients needed
                    !
                    do k=1, nt
                        do n=1, nlat
                            do m=1, mmax
                                brlap(m, n, k) = ZERO
                                bilap(m, n, k) = ZERO
                                crlap(m, n, k) = ZERO
                                cilap(m, n, k) = ZERO
                            end do
                        end do
                        do n=2, nlat
                            brlap(1, n, k) = fnn(n)*br(1, n, k)
                            bilap(1, n, k) = fnn(n)*bi(1, n, k)
                            crlap(1, n, k) = fnn(n)*cr(1, n, k)
                            cilap(1, n, k) = fnn(n)*ci(1, n, k)
                        end do
                        do m=2, mmax
                            do n=m, nlat
                                brlap(m, n, k) = fnn(n)*br(m, n, k)
                                bilap(m, n, k) = fnn(n)*bi(m, n, k)
                                crlap(m, n, k) = fnn(n)*cr(m, n, k)
                                cilap(m, n, k) = fnn(n)*ci(m, n, k)
                            end do
                        end do
                    end do
                case (1, 4, 7)
                    !
                    !     vorticity is zero so cr, ci=0 not used
                    !
                    do k=1, nt
                        do n=1, nlat
                            do m=1, mmax
                                brlap(m, n, k) = ZERO
                                bilap(m, n, k) = ZERO
                            end do
                        end do
                        do n=2, nlat
                            brlap(1, n, k) = fnn(n)*br(1, n, k)
                            bilap(1, n, k) = fnn(n)*bi(1, n, k)
                        end do
                        do m=2, mmax
                            do n=m, nlat
                                brlap(m, n, k) = fnn(n)*br(m, n, k)
                                bilap(m, n, k) = fnn(n)*bi(m, n, k)
                            end do
                        end do
                    end do
                case default
                    !
                    !     divergence is zero so br, bi=0 not used
                    !
                    do k=1, nt
                        do n=1, nlat
                            do m=1, mmax
                                crlap(m, n, k) = ZERO
                                cilap(m, n, k) = ZERO
                            end do
                        end do
                        do n=2, nlat
                            crlap(1, n, k) = fnn(n)*cr(1, n, k)
                            cilap(1, n, k) = fnn(n)*ci(1, n, k)
                        end do
                        do m=2, mmax
                            do n=m, nlat
                                crlap(m, n, k) = fnn(n)*cr(m, n, k)
                                cilap(m, n, k) = fnn(n)*ci(m, n, k)
                            end do
                        end do
                    end do
            end select
        end associate

    end subroutine perform_setup_for_vector_laplacian

    pure function get_workspace_indices_for_inversion(ityp, nlat, mn, lwork) &
        result (return_values)

        ! Dummy arguments
        integer(ip), intent(in) :: ityp
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: mn
        integer(ip), intent(in) :: lwork
        integer(ip)             :: return_values(7)

        associate (i => return_values)
            i(1) = 1
            select case(ityp)
                case(0, 3, 6)
                    i(2) = i(1) + mn
                    i(3) = i(2) + mn
                    i(4) = i(3) + mn
                case(1, 4, 7)
                    i(2) = i(1) + mn
                    i(3) = i(2) + mn
                    i(4) = i(3)
                case default
                    i(2) = i(1)
                    i(3) = i(2) + mn
                    i(4) = i(3) + mn
            end select

            i(5) = i(4) + mn
            i(6) = i(5) + nlat

            select case(ityp)
                case(0, 3, 6)
                    i(7) = lwork - (4 * mn) - nlat
                case default
                    i(7) = lwork - (2 * mn) - nlat
            end select
        end associate

    end function get_workspace_indices_for_inversion

    pure subroutine perform_setup_for_inversion( &
        ityp,  br, bi, cr, ci, brvw, bivw, crvw, civw, fnn)

        ! Dummy arguments
        integer(ip), intent(out) :: ityp
        real(wp),    intent(in)  :: br(:, :, :)
        real(wp),    intent(in)  :: bi(:, :, :)
        real(wp),    intent(in)  :: cr(:, :, :)
        real(wp),    intent(in)  :: ci(:, :, :)
        real(wp),    intent(out) :: brvw(:, :, :)
        real(wp),    intent(out) :: bivw(:, :, :)
        real(wp),    intent(out) :: crvw(:, :, :)
        real(wp),    intent(out) :: civw(:, :, :)
        real(wp),    intent(out) :: fnn(:)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            mmax => size(brvw, dim=1), &
            nlat => size(brvw, dim=2), &
            nt => size(brvw, dim=3) &
           )

            ! Preset coefficient multiplyers in vector
            call compute_coefficient_multipliers(fnn)

            ! Set (u, v) coefficients from br, bi, cr, ci
            select case (ityp)
                case (0, 3, 6)
                    ! All coefficients needed
                    do k=1, nt
                        do n=1, nlat
                            do m=1, mmax
                                brvw(m, n, k) = ZERO
                                bivw(m, n, k) = ZERO
                                crvw(m, n, k) = ZERO
                                civw(m, n, k) = ZERO
                            end do
                        end do
                        do n=2, nlat
                            brvw(1, n, k) = br(1, n, k)/fnn(n)
                            bivw(1, n, k) = bi(1, n, k)/fnn(n)
                            crvw(1, n, k) = cr(1, n, k)/fnn(n)
                            civw(1, n, k) = ci(1, n, k)/fnn(n)
                        end do
                        do m=2, mmax
                            do n=m, nlat
                                brvw(m, n, k) = br(m, n, k)/fnn(n)
                                bivw(m, n, k) = bi(m, n, k)/fnn(n)
                                crvw(m, n, k) = cr(m, n, k)/fnn(n)
                                civw(m, n, k) = ci(m, n, k)/fnn(n)
                            end do
                        end do
                    end do
                case (1, 4, 7)
                    ! Vorticity is zero so cr, ci=0 not used
                    do k=1, nt
                        do n=1, nlat
                            do m=1, mmax
                                brvw(m, n, k) = ZERO
                                bivw(m, n, k) = ZERO
                            end do
                        end do
                        do n=2, nlat
                            brvw(1, n, k) = br(1, n, k)/fnn(n)
                            bivw(1, n, k) = bi(1, n, k)/fnn(n)
                        end do
                        do m=2, mmax
                            do n=m, nlat
                                brvw(m, n, k) = br(m, n, k)/fnn(n)
                                bivw(m, n, k) = bi(m, n, k)/fnn(n)
                            end do
                        end do
                    end do
                case default
                    ! Divergence is zero so br, bi=0 not used
                    do k=1, nt
                        do n=1, nlat
                            do m=1, mmax
                                crvw(m, n, k) = ZERO
                                civw(m, n, k) = ZERO
                            end do
                        end do
                        do n=2, nlat
                            crvw(1, n, k) = cr(1, n, k)/fnn(n)
                            civw(1, n, k) = ci(1, n, k)/fnn(n)
                        end do
                        do m=2, mmax
                            do n=m, nlat
                                crvw(m, n, k) = cr(m, n, k)/fnn(n)
                                civw(m, n, k) = ci(m, n, k)/fnn(n)
                            end do
                        end do
                    end do
            end select
        end associate

    end subroutine perform_setup_for_inversion

end module vector_laplacian_routines
