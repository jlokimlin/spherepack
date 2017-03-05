module divergence_routines

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
    public :: divec, dives, divgc, divgs
    public :: idivec, idives, idivgc, idivgs
    public :: perform_setup_for_inversion
    public :: perform_setup_for_divergence

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: SQRT2 = sqrt(TWO)

    ! Declare interfaces for submodule implementation
    interface
        module subroutine divgs(nlat, nlon, isym, nt, divg, idiv, jdiv, br, bi, mdb, ndb, &
            wshsgs, lshsgs, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: divg(idiv, jdiv, nt)
            integer(ip), intent(in)  :: idiv
            integer(ip), intent(in)  :: jdiv
            real(wp),    intent(in)  :: br(mdb, ndb, nt)
            real(wp),    intent(in)  :: bi(mdb, ndb, nt)
            integer(ip), intent(in)  :: mdb
            integer(ip), intent(in)  :: ndb
            real(wp),    intent(in)  :: wshsgs(lshsgs)
            integer(ip), intent(in)  :: lshsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine divgs

        module subroutine divgc(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            wshsgc, lshsgc, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: dv(idv, jdv, nt)
            integer(ip), intent(in)  :: idv
            integer(ip), intent(in)  :: jdv
            real(wp),    intent(in)  :: br(mdb, ndb, nt)
            real(wp),    intent(in)  :: bi(mdb, ndb, nt)
            integer(ip), intent(in)  :: mdb
            integer(ip), intent(in)  :: ndb
            real(wp),    intent(in)  :: wshsgc(lshsgc)
            integer(ip), intent(in)  :: lshsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine divgc

        module subroutine dives(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            wshses, lshses, work, lwork, ierror)
            
            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: dv(idv, jdv, nt)
            integer(ip), intent(in)  :: idv
            integer(ip), intent(in)  :: jdv
            real(wp),    intent(in)  :: br(mdb, ndb, nt)
            real(wp),    intent(in)  :: bi(mdb, ndb, nt)
            integer(ip), intent(in)  :: mdb
            integer(ip), intent(in)  :: ndb
            real(wp),    intent(in)  :: wshses(lshses)
            integer(ip), intent(in)  :: lshses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine dives

        module subroutine divec(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            wshsec, lshsec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: dv(idv, jdv, nt)
            integer(ip), intent(in)  :: idv
            integer(ip), intent(in)  :: jdv
            real(wp),    intent(in)  :: br(mdb, ndb, nt)
            real(wp),    intent(in)  :: bi(mdb, ndb, nt)
            integer(ip), intent(in)  :: mdb
            integer(ip), intent(in)  :: ndb
            real(wp),    intent(in)  :: wshsec(lshsec)
            integer(ip), intent(in)  :: lshsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine divec

        module subroutine idivec(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
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
        end subroutine idivec

        module subroutine idives(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
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
        end subroutine idives

        module subroutine idivgc(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
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
        end subroutine idivgc

        module subroutine idivgs(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
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
        end subroutine idivgs
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

    pure subroutine perform_setup_for_divergence(nlon, a, b, br, bi, sqnn)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: a(:, :, :)
        real(wp),    intent(out) :: b(:, :, :)
        real(wp),    intent(in)  :: br(:, :, :)
        real(wp),    intent(in)  :: bi(:, :, :)
        real(wp),    intent(out) :: sqnn(:)

        ! Local variables
        integer(ip) :: k, n, m, mmax

        associate (&
            nlat => size(sqnn), &
            nt => size(a, dim=3) &
           )

            ! Set coefficient multiplyers
            call compute_coefficient_multipliers(sqnn)

            ! Compute divergence scalar coefficients for each vector field
            do k=1, nt
                a(:, :, k) = ZERO
                b(:, :, k) = ZERO

                ! Compute m=0 coefficients
                do n=2, nlat
                    a(1, n, k) = -sqnn(n)*br(1, n, k)
                    b(1, n, k) = -sqnn(n)*bi(1, n, k)
                end do

                ! Compute m > 0 coefficients using vector spherepack value for mmax
                mmax = min(nlat, (nlon+1)/2)
                do m=2, mmax
                    do n=m, nlat
                        a(m, n, k) = -sqnn(n)*br(m, n, k)
                        b(m, n, k) = -sqnn(n)*bi(m, n, k)
                    end do
                end do
            end do
        end associate

    end subroutine perform_setup_for_divergence

    pure subroutine perform_setup_for_inversion(isym, ityp, a, b, sqnn, pertrb, br, bi)

        ! Dummy arguments
        integer(ip), intent(in)  :: isym
        integer(ip), intent(out) :: ityp
        real(wp),    intent(in)  :: a(:, :, :)
        real(wp),    intent(in)  :: b(:, :, :)
        real(wp),    intent(out) :: sqnn(:)
        real(wp),    intent(out) :: pertrb(:)
        real(wp),    intent(out) :: br(:, :, :)
        real(wp),    intent(out) :: bi(:, :, :)

        ! Local variables
        integer(ip) :: k, n, m

        associate (&
            mmax => size(br, dim=1), &
            nlat => size(br, dim=2), &
            nt => size(br, dim=3) &
           )

            ! Preset coefficient multiplyers in vector
            call compute_coefficient_multipliers(sqnn)

            ! Compute multiple vector fields coefficients
            do k=1, nt

                ! Set divergence field perturbation adjustment
                pertrb(k) = get_perturbation(a, k)

                ! Preset br, bi to 0.0
                do n=1, nlat
                    do m=1, mmax
                        br(m, n, k) = ZERO
                        bi(m, n, k) = ZERO
                    end do
                end do

                ! Compute m=0 coefficients
                do n=2, nlat
                    br(1, n, k) = -a(1, n, k)/sqnn(n)
                    bi(1, n, k) = -b(1, n, k)/sqnn(n)
                end do

                ! Compute m > 0 coefficients
                do m=2, mmax
                    do n=m, nlat
                        br(m, n, k) = -a(m, n, k)/sqnn(n)
                        bi(m, n, k) = -b(m, n, k)/sqnn(n)
                    end do
                end do
            end do

            ! Set ityp for vector synthesis with curl=0
            select case (isym)
                case (0)
                    ityp = 1
                case (1)
                    ityp = 4
                case (2)
                    ityp = 7
            end select
        end associate

    end subroutine perform_setup_for_inversion

end module divergence_routines
