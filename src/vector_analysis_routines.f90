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
            mdab, ndab, wvhagc, lvhagc, work, lwork, ierror)

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
            real(wp),    intent(in)  :: wvhagc(lvhagc)
            integer(ip), intent(in)  :: lvhagc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine vhagc

        module subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mdab, ndab, wvhags, lvhags, work, lwork, ierror)

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
            real(wp),    intent(in)  :: wvhags(lvhags)
            integer(ip), intent(in)  :: lvhags
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
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
        procedure, nopass :: get_lwork => get_vhaes_lwork
        procedure, nopass :: get_ldwork => get_vhaes_ldwork
        procedure, nopass :: get_legendre_workspace_size => &
            get_vhaes_legendre_workspace_size
    end type VhaesAux

    type, public :: VhagsAux
    contains
        ! Type-bound procedures
        procedure, nopass :: vhags
        procedure, nopass :: vhagsi
        procedure, nopass :: get_lvhags
        procedure, nopass :: get_ldwork
        procedure, nopass :: get_legendre_workspace_size
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

    pure function get_ldwork(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value


        return_value = (3*nlat*(nlat+3)+2)/2

    end function get_ldwork

    pure function get_legendre_workspace_size(nlat, nlon, nt, ityp) result (return_value)

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
        if (ityp <= 2) then
            ! Set workspace size
            return_value = max(3*nlat*(nlat+1)+2, (2*nt_op+1)*nlat*nlon)
        else
            ! Compute parity
            select case (mod(nlat, 2))
                case (0)
                    l2 = nlat/2
                case default
                    l2 = (nlat + 1)/2
            end select
            ! Set workspace size
            return_value = max(3*nlat*(nlat+1)+2, (2*nt_op+1)*l2*nlon)
        end if

    end function get_legendre_workspace_size

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

    pure function get_vhaes_lwork(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackUtility) :: util


        call util%compute_parity(nlat, nlon, l1, l2)

        return_value = 3*(max(l1-2, 0)*(2*nlat-l1-1))/2+5*l2*nlat

    end function get_vhaes_lwork

    pure function get_vhaes_ldwork(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value


        return_value = 2*(nlat+1)

    end function get_vhaes_ldwork

    pure function get_vhaes_legendre_workspace_size(nlat, nlon, nt, ityp) result (return_value)

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

    end function get_vhaes_legendre_workspace_size

end module vector_analysis_routines
