module divergence_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsec, &
        shses, &
        shsgc, &
        shsgs

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: divec, dives, divgc, divgs

    interface
        module subroutine divgs(nlat, nlon, isym, nt, divg, idiv, jdiv, br, bi, mdb, ndb, &
            wshsgs, lshsgs, work, lwork, ierror)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
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
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer(ip), intent(in)     :: nlat
            integer(ip), intent(in)     :: nlon
            integer(ip), intent(in)     :: isym
            integer(ip), intent(in)     :: nt
            real(wp),    intent(inout)  :: dv(idv, jdv,*)
            integer(ip), intent(in)     :: idv
            integer(ip), intent(in)     :: jdv
            real(wp),    intent(inout)  :: br(mdb, ndb,*)
            real(wp),    intent(inout)  :: bi(mdb, ndb,*)
            integer(ip), intent(in)     :: mdb
            integer(ip), intent(in)     :: ndb
            real(wp),    intent(inout)  :: wshsec(lshsec)
            integer(ip), intent(in)     :: lshsec
            real(wp),    intent(inout)  :: work(lwork)
            integer(ip), intent(in)     :: lwork
            integer(ip), intent(out)    :: ierror
        end subroutine divec
    end interface

    !------------------------------------------------------------------
    ! Parameters confined to the module
    !------------------------------------------------------------------
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp
    !------------------------------------------------------------------

end module divergence_routines
