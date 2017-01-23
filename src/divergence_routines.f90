module divergence_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsec, &
        shses, &
        shsgc, &
        shsgs

    use vector_synthesis_routines, only: &
        vhsgc, &
        vhses, &
        vhsec, &
        vhsgs

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: divec, dives, divgc, divgs
    public :: idivec, idives, idivgc, idivgs
    public :: compute_coefficient_multipliers
    public :: get_perturbation

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
            integer(ip), intent(in)     :: nlat
            integer(ip), intent(in)     :: nlon
            integer(ip), intent(in)     :: isym
            integer(ip), intent(in)     :: nt
            real(wp),    intent(inout)  :: divg(idiv, jdiv, nt)
            integer(ip), intent(in)     :: idiv
            integer(ip), intent(in)     :: jdiv
            real(wp),    intent(inout)  :: br(mdb, ndb, nt)
            real(wp),    intent(inout)  :: bi(mdb, ndb, nt)
            integer(ip), intent(in)     :: mdb
            integer(ip), intent(in)     :: ndb
            real(wp),    intent(inout)  :: wshsgs(lshsgs)
            integer(ip), intent(in)     :: lshsgs
            real(wp),    intent(inout)  :: work(lwork)
            integer(ip), intent(in)     :: lwork
            integer(ip), intent(out)    :: ierror
        end subroutine divgs

        module subroutine divgc(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            wshsgc, lshsgc, work, lwork, ierror)
            real(wp) :: dv(idv, jdv, nt), br(mdb, ndb, nt), bi(mdb, ndb, nt)
            integer(ip) :: idv
            integer(ip) :: ierror
            integer(ip) :: isym
            integer(ip) :: jdv
            integer(ip) :: lshsgc
            integer(ip) :: lwork
            integer(ip) :: mab
            integer(ip) :: mdb
            integer(ip) :: ndb
            integer(ip) :: nlat
            integer(ip) :: nlon
            integer(ip) :: nt
            real(wp) :: wshsgc(lshsgc), work(lwork)
        end subroutine divgc

        module subroutine dives(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            wshses, lshses, work, lwork, ierror)
            real(wp) :: dv(idv, jdv, nt), br(mdb, ndb, nt), bi(mdb, ndb, nt)
            integer(ip) :: idv
            integer(ip) :: ierror
            integer(ip) :: isym
            integer(ip) :: iwk
            integer(ip) :: jdv
            integer(ip) :: lshses
            integer(ip) :: lwork
            integer(ip) :: mab
            integer(ip) :: mdb
            integer(ip) :: mmax
            integer(ip) :: ndb
            integer(ip) :: nlat
            integer(ip) :: nlon
            integer(ip) :: nt
            real(wp) :: wshses(lshses), work(lwork)
        end subroutine dives

        module subroutine divec(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            wshsec, lshsec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)     :: nlat
            integer(ip), intent(in)     :: nlon
            integer(ip), intent(in)     :: isym
            integer(ip), intent(in)     :: nt
            real(wp),    intent(inout)  :: dv(idv, jdv, nt)
            integer(ip), intent(in)     :: idv
            integer(ip), intent(in)     :: jdv
            real(wp),    intent(inout)  :: br(mdb, ndb, nt)
            real(wp),    intent(inout)  :: bi(mdb, ndb, nt)
            integer(ip), intent(in)     :: mdb
            integer(ip), intent(in)     :: ndb
            real(wp),    intent(inout)  :: wshsec(lshsec)
            integer(ip), intent(in)     :: lshsec
            real(wp),    intent(inout)  :: work(lwork)
            integer(ip), intent(in)     :: lwork
            integer(ip), intent(out)    :: ierror
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

        associate( nlat => size(sqnn) )
            do n=2, nlat
                fn = real(n - 1, kind=wp)
                sqnn(n) = sqrt(fn * (fn + ONE))
            end do
        end associate

    end subroutine compute_coefficient_multipliers

    pure function get_perturbation(a, k) &
        result(return_value)

        ! Dummy arguments
        real(wp),    intent(in) :: a(:,:,:)
        integer(ip), intent(in) :: k
        real(wp)                :: return_value

        return_value = a(1,1,k)/(TWO * SQRT2)

    end function get_perturbation

end module divergence_routines
