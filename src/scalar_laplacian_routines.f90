module scalar_laplacian_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsec, shses, shsgc, shsgs

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: slapec, slapes, slapgc, slapgs
    public :: islapec, islapes, islapgc, islapgs
    public :: perform_setup_for_inversion
    public :: perform_setup_for_scalar_laplacian

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine slapec(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshsec, lshsec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsec(lshsec)
            integer(ip), intent(in)  :: lshsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine slapec

        module subroutine slapes(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshses, lshses, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshses(lshses)
            integer(ip), intent(in)  :: lshses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine slapes

        module subroutine slapgc(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshsgc, lshsgc, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgc(lshsgc)
            integer(ip), intent(in)  :: lshsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine slapgc

        module subroutine slapgs(nlat, nlon, isym, nt, slap, ids, jds, a, b, mdab, ndab, &
            wshsgs, lshsgs, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: slap(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgs(lshsgs)
            integer(ip), intent(in)  :: lshsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine slapgs

        module subroutine islapec(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshsec, lshsec, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(nt)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsec(lshsec)
            integer(ip), intent(in)  :: lshsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine islapec

        module subroutine islapes(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshses, lshses, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(nt)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshses(lshses)
            integer(ip), intent(in)  :: lshses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine islapes

        module subroutine islapgc(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshsgc, lshsgc, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(nt)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgc(lshsgc)
            integer(ip), intent(in)  :: lshsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine islapgc

        module subroutine islapgs(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
            mdab, ndab, wshsgs, lshsgs, work, lwork, pertrb, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(in)  :: xlmbda(nt)
            real(wp),    intent(out) :: sf(ids, jds, nt)
            integer(ip), intent(in)  :: ids
            integer(ip), intent(in)  :: jds
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgs(lshsgs)
            integer(ip), intent(in)  :: lshsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: pertrb(nt)
            integer(ip), intent(out) :: ierror
        end subroutine islapgs
    end interface

contains

    pure subroutine compute_coefficient_multipliers(fnn)

        ! Dummy arguments
        real(wp), intent(out) :: fnn(:)

        ! Local variables
        integer(ip) :: n
        real(wp)    :: fn

        associate( nlat => size(fnn) )
            do n=2, nlat
                fn = real(n - 1, kind=wp)
                fnn(n) = fn * (fn + ONE)
            end do
        end associate

    end subroutine compute_coefficient_multipliers

    pure subroutine perform_setup_for_scalar_laplacian(a, b, alap, blap, fnn)

        ! Dummy arguments
        real(wp),    intent(in)  :: a(:, :, :)
        real(wp),    intent(in)  :: b(:, :, :)
        real(wp),    intent(out) :: alap(:, :, :)
        real(wp),    intent(out) :: blap(:, :, :)
        real(wp),    intent(out) :: fnn(:)

        ! Local variables
        integer(ip) :: k, n, m

        associate( &
            mmax => size(alap, dim=1), &
            nlat => size(alap, dim=2), &
            nt => size(alap, dim=3) &
            )

            ! Set coefficient multiplyers
            call compute_coefficient_multipliers(fnn)

            ! Compute scalar laplacian coefficients for each vector field
            do k=1, nt
                alap(:, :, k) = ZERO
                blap(:, :, k) = ZERO

                ! Compute m = 0 coefficients
                do n=2, nlat
                    alap(1, n, k) = -fnn(n) * a(1, n, k)
                    blap(1, n, k) = -fnn(n) * b(1, n, k)
                end do

                ! Compute m > 0 coefficients
                do m=2, mmax
                    do n=m, nlat
                        alap(m, n, k) = -fnn(n)*a(m, n, k)
                        blap(m, n, k) = -fnn(n)*b(m, n, k)
                    end do
                end do
            end do
        end associate

    end subroutine perform_setup_for_scalar_laplacian

    pure subroutine perform_setup_for_inversion(a, b, as, bs, fnn, xlmbda, pertrb)

        ! Dummy arguments
        real(wp), intent(in)  :: a(:, :, :)
        real(wp), intent(in)  :: b(:, :, :)
        real(wp), intent(out) :: as(:, :, :)
        real(wp), intent(out) :: bs(:, :, :)
        real(wp), intent(out) :: fnn(:)
        real(wp), intent(in)  :: xlmbda(:)
        real(wp), intent(out) :: pertrb(:)

        ! Local variables
        integer(ip) :: k, n, m

        associate( &
            mmax => size(as, dim=1), &
            nlat => size(as, dim=2), &
            nt => size(as, dim=3) &
            )

            ! Preset coefficient multiplyers
            call compute_coefficient_multipliers(fnn)

            ! Preset synthesis coefficients to zero
            as = ZERO
            bs = ZERO

            do k=1, nt

                ! Compute synthesis coefficients for xlmbda zero or nonzero
                if (xlmbda(k) == ZERO) then
                    do n=2, nlat
                        as(1, n, k) = -a(1, n, k)/fnn(n)
                        bs(1, n, k) = -b(1, n, k)/fnn(n)
                    end do
                    do m=2, mmax
                        do n=m, nlat
                            as(m, n, k) = -a(m, n, k)/fnn(n)
                            bs(m, n, k) = -b(m, n, k)/fnn(n)
                        end do
                    end do
                else
                    ! xlmbda nonzero so operator invertible unless
                    ! -n*(n-1) = xlmbda(k) < 0.0  for some n
                    !
                    pertrb(k) = ZERO

                    do n=1, nlat
                        as(1, n, k) = -a(1, n, k)/(fnn(n) + xlmbda(k))
                        bs(1, n, k) = -b(1, n, k)/(fnn(n) + xlmbda(k))
                    end do

                    do m=2, mmax
                        do n=m, nlat
                            as(m, n, k) = -a(m, n, k)/(fnn(n) + xlmbda(k))
                            bs(m, n, k) = -b(m, n, k)/(fnn(n) + xlmbda(k))
                        end do
                    end do
                end if
            end do
        end associate

    end subroutine perform_setup_for_inversion

end module scalar_laplacian_routines
