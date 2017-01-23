module vector_synthesis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_SpherepackAux, only: &
        SpherepackAux

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vhsgc, vhsgci
    public :: vhses, vhsesi, VhsesAux
    public :: vhsec, vhseci
    public :: vhsgs, vhsgsi, VhsgsAux
    
    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    ! Declare interfaces for submodule implementation
    interface
        module subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhses, lvhses, work, lwork, ierror)

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
            real(wp),    intent(in)  :: wvhses(lvhses)
            integer(ip), intent(in)  :: lvhses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror

        end subroutine vhses

        module subroutine vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, &
            ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhses(lvhses)
            integer(ip), intent(in)  :: lvhses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror

        end subroutine vhsesi

        module subroutine vhsgs(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
            mdab,ndab,wvhsgs,lvhsgs,work,lwork,ierror)

            ! Dummy arguments
            integer(ip), intent(in)     :: nlat
            integer(ip), intent(in)     :: nlon
            integer(ip), intent(in)     :: ityp
            integer(ip), intent(in)     :: nt
            real(wp),    intent(out)    :: v(idvw,jdvw,nt)
            real(wp),    intent(out)    :: w(idvw,jdvw,nt)
            integer(ip), intent(in)     :: idvw
            integer(ip), intent(in)     :: jdvw
            real(wp),    intent(in)     :: br(mdab,ndab,nt)
            real(wp),    intent(in)     :: bi(mdab,ndab,nt)
            real(wp),    intent(in)     :: cr(mdab,ndab,nt)
            real(wp),    intent(in)     :: ci(mdab,ndab,nt)
            integer(ip), intent(in)     :: mdab
            integer(ip), intent(in)     :: ndab
            real(wp),    intent(in)     :: wvhsgs(lvhsgs)
            integer(ip), intent(in)     :: lvhsgs
            real(wp),    intent(out)    :: work(lwork)
            integer(ip), intent(in)     :: lwork
            integer(ip), intent(out)    :: ierror

        end subroutine vhsgs

        module subroutine vhsgsi(nlat,nlon,wvhsgs,lvhsgs,dwork,ldwork,ierror)

            ! Dummy arguments
            integer(ip), intent(in)     :: nlat
            integer(ip), intent(in)     :: nlon
            real(wp),    intent(out)    :: wvhsgs(lvhsgs)
            integer(ip), intent(in)     :: lvhsgs
            real(wp),    intent(out)    :: dwork(ldwork)
            integer(ip), intent(in)     :: ldwork
            integer(ip), intent(out)    :: ierror

        end subroutine vhsgsi

        module subroutine vhsgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsgc, lvhsgc, work, lwork, ierror)

            ! Dummy arguments
            real(wp) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
            real(wp) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)
            integer(ip) :: idvw
            integer(ip) :: ierror
            integer(ip) :: ityp
            integer(ip) :: jdvw
            integer(ip) :: lvhsgc
            integer(ip) :: lwork
            integer(ip) :: mdab
            integer(ip) :: ndab
            integer(ip) :: nlat
            integer(ip) :: nlon
            integer(ip) :: nt
            real(wp) :: v(idvw, jdvw, nt), w(idvw, jdvw, nt)
            real(wp) :: work(lwork), wvhsgc(lvhsgc)
        end subroutine vhsgc

        module subroutine vhsgci(nlat, nlon, wvhsgc, lvhsgc, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip) :: ierror
            integer(ip) :: ldwork
            integer(ip) :: lvhsgc
            integer(ip) :: nlat
            integer(ip) :: nlon
            real(wp)    :: wvhsgc(lvhsgc)
            real(wp)    :: dwork(ldwork)
        end subroutine vhsgci

        module subroutine vhsec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhsec, lvhsec, work, lwork, ierror)

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
            real(wp),    intent(in)  :: wvhsec(lvhsec)
            integer(ip), intent(in)  :: lvhsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vhsec

        module subroutine vhseci(nlat, nlon, wvhsec, lvhsec, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wvhsec(lvhsec)
            integer(ip), intent(in)  :: lvhsec
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine vhseci
    end interface

    type, public :: VhsesAux
    contains
        ! Type-bound procedures
        procedure, nopass :: vhses
        procedure, nopass :: vhsesi
        procedure, nopass :: get_lvhses
        procedure, nopass :: get_lwork => get_vhses_lwork
        procedure, nopass :: get_ldwork => get_vhses_ldwork
        procedure, nopass :: get_legendre_workspace_size &
            => get_vhses_legendre_workspace_size
    end type VhsesAux

    type, public :: VhsgsAux
    contains
        ! Type-bound procedures
        procedure, nopass :: vhsgs
        procedure, nopass :: vhsgsi
        procedure, nopass :: get_lvhsgs
        procedure, nopass :: get_ldwork
        procedure, nopass :: get_legendre_workspace_size
    end type VhsgsAux

contains

    pure function get_lvhsgs(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables
        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = max(l1*l2*(2*nlat-l1+1)+nlon+15+2*nlat, 5*(nlat**2)*nlon)

    end function get_lvhsgs

    pure function get_ldwork(nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value

        return_value = (3*nlat*(nlat+3)+2)/2

    end function get_ldwork

    pure function get_legendre_workspace_size(nlat, nlon, nt, ityp) &
        result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat
        integer(ip),           intent(in) :: nlon
        integer(ip), optional, intent(in) :: ityp
        integer(ip), optional, intent(in) :: nt
        integer(ip)                        :: return_value

        ! Local variables
        integer(ip) :: nt_op, ityp_op, l2

        !  Address optional arguments
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        if (present(ityp)) then
            ityp_op = ityp
        else
            ityp_op = 1
        end if

        !
        !  Compute workspace size
        !
        if (ityp <= 2) then
            ! Set workspace size
            return_value = (2*nt_op+1)*nlat*nlon
        else
            ! Compute parity
            select case (mod(nlat, 2))
                case (0)
                    l2 = nlat/2
                case default
                    l2 = (nlat + 1)/2
            end select
            ! Set workspace size
            return_value = (2*nt+1)*l2*nlon
        end if

    end function get_legendre_workspace_size

    pure function get_lvhses(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux


        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = l1*l2*(2*nlat-l1+1)+nlon+15

    end function get_lvhses

    pure function get_vhses_lwork(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux


        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = 3*(max(l1-2, 0)*(2*nlat-l1-1))/2+5*l2*nlat

    end function get_vhses_lwork

    pure function get_vhses_ldwork(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value


        return_value = 2*(nlat+1)

    end function get_vhses_ldwork

    pure function get_vhses_legendre_workspace_size(nlat, nlon, nt, ityp) result (return_value)

        ! Dummy arguments

        integer(ip),           intent(in) :: nlat
        integer(ip),           intent(in) :: nlon
        integer(ip), optional, intent(in) :: nt
        integer(ip), optional, intent(in) :: ityp
        integer(ip)                        :: return_value

        ! Local variables

        integer(ip) :: nt_op, ityp_op, l2


        !
        !  Address optional arguments
        !
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        if (present(ityp)) then
            ityp_op = ityp
        else
            ityp_op = 0
        end if

        !
        !  Compute workspace size
        !
        if (ityp_op <= 2) then
            ! Set workspace size
            return_value = (2*nt_op+1)*nlat*nlon
        else
            ! Compute parity
            select case (mod(nlat, 2))
                case (0)
                    l2 = nlat/2
                case default
                    l2 = (nlat + 1)/2
            end select
            ! Set workspace size
            return_value = (2*nt_op+1)*l2*nlon
        end if

    end function get_vhses_legendre_workspace_size

end module vector_synthesis_routines
