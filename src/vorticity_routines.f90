module vorticity_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsec, shses, shsgc, shsgs

    use vector_synthesis_routines, only: &
        vhses, vhsec, vhsgc, vhsgs

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vrtec, vrtes, vrtgc, vrtgs
    public :: ivrtec, ivrtes, ivrtgc, ivrtgs
    public :: perform_setup_for_inversion
    public :: perform_setup_for_vorticity

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: SQRT2 = sqrt(TWO)

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vrtec(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshsec, lshsec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshsec(lshsec)
            integer(ip), intent(in)  :: lshsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vrtec

        module subroutine vrtes(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshses, lshses, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshses(lshses)
            integer(ip), intent(in)  :: lshses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vrtes

        module subroutine vrtgc(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshsgc, lshsgc, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshsgc(lshsgc)
            integer(ip), intent(in)  :: lshsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vrtgc

        module subroutine vrtgs(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            wshsgs, lshsgs, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
            integer(ip), intent(in)  :: ivrt
            integer(ip), intent(in)  :: jvrt
            real(wp),    intent(in)  :: cr(mdc, ndc, nt)
            real(wp),    intent(in)  :: ci(mdc, ndc, nt)
            integer(ip), intent(in)  :: mdc
            integer(ip), intent(in)  :: ndc
            real(wp),    intent(in)  :: wshsgs(lshsgs)
            integer(ip), intent(in)  :: lshsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vrtgs

        module subroutine ivrtec(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhsec, lvhsec, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhsec(lvhsec)
            integer(ip), intent(in)  :: lvhsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtec

        module subroutine ivrtes(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhses, lvhses, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhses(lvhses)
            integer(ip), intent(in)  :: lvhses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtes

        module subroutine ivrtgc(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhsgc, lvhsgc, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhsgc(lvhsgc)
            integer(ip), intent(in)  :: lvhsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtgc

        module subroutine ivrtgs(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
            wvhsgs, lvhsgs, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: v(idvw, jdvw, nt)
            real(wp),    intent(out) :: w(idvw, jdvw, nt)
            integer(ip), intent(in)  :: idvw
            integer(ip), intent(in)  :: jdvw
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(out) :: wvhsgs(lvhsgs)
            integer(ip), intent(in)  :: lvhsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine ivrtgs
    end interface

contains

    pure subroutine compute_coefficient_multipliers(sqnn)

        ! Dummy arguments
        real(wp), intent(out) :: sqnn(:)

        ! Local variables
        integer(ip) :: n
        real(wp)    :: fn

        associate (nlat => size(sqnn))
            do n=2, nlat
                fn = real(n - 1, kind=wp)
                sqnn(n) = sqrt(fn * (fn + ONE))
            end do
        end associate

    end subroutine compute_coefficient_multipliers

    pure function get_perturbation(a, k) &
        result(return_value)

        ! Dummy arguments
        real(wp),    intent(in) :: a(:, :, :)
        integer(ip), intent(in) :: k
        real(wp)                :: return_value

        return_value = a(1, 1, k)/(TWO * SQRT2)

    end function get_perturbation

    pure subroutine perform_setup_for_vorticity(nlon, a, b, cr, ci, sqnn)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: a(:, :, :)
        real(wp),    intent(out) :: b(:, :, :)
        real(wp),    intent(in)  :: cr(:, :, :)
        real(wp),    intent(in)  :: ci(:, :, :)
        real(wp),    intent(out) :: sqnn(:)

        ! Local variables
        integer(ip) :: k, n, m, mmax

        associate (&
            nlat => size(sqnn), &
            nt => size(a, dim=3) &
           )

            ! Set coefficient multiplyers
            call compute_coefficient_multipliers(sqnn)

            ! Compute vorticity scalar coefficients for each vector field
            do k=1, nt
                a(:, :, k) = ZERO
                b(:, :, k) = ZERO

                ! Compute m=0 coefficients
                do n=2, nlat
                    a(1, n, k) = sqnn(n) * cr(1, n, k)
                    b(1, n, k) = sqnn(n) * ci(1, n, k)
                end do

                ! Compute m > 0 coefficients using vector spherepack value for mmax
                mmax = min(nlat, (nlon+1)/2)
                do m=2, mmax
                    do n=m, nlat
                        a(m, n, k) = sqnn(n) * cr(m, n, k)
                        b(m, n, k) = sqnn(n) * ci(m, n, k)
                    end do
                end do
            end do
        end associate

    end subroutine perform_setup_for_vorticity

    pure subroutine perform_setup_for_inversion(isym, ityp, a, b, sqnn, pertrb, cr, ci)

        ! Dummy arguments
        integer(ip), intent(in)  :: isym
        integer(ip), intent(out) :: ityp
        real(wp),    intent(in)  :: a(:, :, :)
        real(wp),    intent(in)  :: b(:, :, :)
        real(wp),    intent(out) :: sqnn(:)
        real(wp),    intent(out) :: pertrb(:)
        real(wp),    intent(out) :: cr(:, :, :)
        real(wp),    intent(out) :: ci(:, :, :)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            mmax => size(cr, dim=1), &
            nlat => size(cr, dim=2), &
            nt => size(cr, dim=3) &
           )

            ! Preset coefficient multiplyers in vector
            call compute_coefficient_multipliers(sqnn)

            ! Compute multiple vector fields coefficients
            do k=1, nt

                ! Set vorticity field perturbation adjustment
                pertrb(k) = get_perturbation(a, k)

                ! Preset cr, ci to 0.0
                cr(:, :, k) = ZERO
                ci(:, :, k) = ZERO

                ! Compute m = 0 coefficients
                do n=2, nlat
                    cr(1, n, k) = a(1, n, k)/sqnn(n)
                    ci(1, n, k) = b(1, n, k)/sqnn(n)
                end do

                ! Compute m > 0 coefficients
                do m=2, mmax
                    do n=m, nlat
                        cr(m, n, k) = a(m, n, k)/sqnn(n)
                        ci(m, n, k) = b(m, n, k)/sqnn(n)
                    end do
                end do
            end do

            ! Set ityp for vector synthesis with divergence=0
            select case (isym)
                case (0)
                    ityp = 2
                case (1)
                    ityp = 5
                case (2)
                    ityp = 8
            end select
        end associate

    end subroutine perform_setup_for_inversion

end module vorticity_routines
