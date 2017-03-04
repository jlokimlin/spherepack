module scalar_synthesis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SpherepackUtility, only: &
        SpherepackUtility

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: shses, shsesi
    public :: shsec, shseci
    public :: shsgs, shsgsi
    public :: shsgc, shsgci
    public :: ShsesAux, ShsgsAux

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp

    type, public :: ShsesAux
    contains
        ! Type-bound procedures
        procedure, nopass :: shses
        procedure, nopass :: shsesi
        procedure, nopass :: get_lshses
        procedure, nopass :: get_lwork &
            => get_lwork_shses
        procedure, nopass :: get_ldwork &
            => get_ldwork_shses
        procedure, nopass :: get_legendre_workspace_size &
            => get_legendre_workspace_size_shses
    end type ShsesAux

    type, public :: ShsgsAux
    contains
        ! Type-bound procedures
        procedure, nopass :: shsgs
        procedure, nopass :: shsgsi
        procedure, nopass :: get_lshsgs
        procedure, nopass :: get_lwork &
            => get_lwork_shsgs
        procedure, nopass :: get_ldwork &
            => get_ldwork_shsgs
        procedure, nopass :: get_legendre_workspace_size &
            => get_legendre_workspace_size_shsgs
    end type ShsgsAux

    ! Declare interfaces for submodule implementation
    interface
        module subroutine shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshses, lshses, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshses(lshses)
            integer(ip), intent(in)  :: lshses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shses

        module subroutine shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, &
            ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshses(lshses)
            integer(ip), intent(in)  :: lshses
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine shsesi

        module subroutine shsgs(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshsgs, lshsgs, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: mode
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgs(lshsgs)
            integer(ip), intent(in)  :: lshsgs
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shsgs

        module subroutine shsgsi(nlat, nlon, wshsgs, lshsgs, work, lwork, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            real(wp),    intent(out)  :: wshsgs(lshsgs)
            integer(ip), intent(in)   :: lshsgs
            real(wp),    intent(out)  :: work(lwork)
            integer(ip), intent(in)   :: lwork
            real(wp),    intent(out)  :: dwork(ldwork)
            integer(ip), intent(in)   :: ldwork
            integer(ip), intent(out)  :: ierror
        end subroutine shsgsi

        module subroutine shsec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshsec, lshsec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: isym
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsec(lshsec)
            integer(ip), intent(in)  :: lshsec
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shsec

        module subroutine shseci(nlat, nlon, wshsec, lshsec, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshsec(lshsec)
            integer(ip), intent(in)  :: lshsec
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine shseci

        module subroutine shsgc(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshsgc, lshsgc, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            integer(ip), intent(in)  :: mode
            integer(ip), intent(in)  :: nt
            real(wp),    intent(out) :: g(idg, jdg, nt)
            integer(ip), intent(in)  :: idg
            integer(ip), intent(in)  :: jdg
            real(wp),    intent(in)  :: a(mdab, ndab, nt)
            real(wp),    intent(in)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)  :: mdab
            integer(ip), intent(in)  :: ndab
            real(wp),    intent(in)  :: wshsgc(lshsgc)
            integer(ip), intent(in)  :: lshsgc
            real(wp),    intent(out) :: work(lwork)
            integer(ip), intent(in)  :: lwork
            integer(ip), intent(out) :: ierror
        end subroutine shsgc

        module subroutine shsgci(nlat, nlon, wshsgc, lshsgc, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshsgc(lshsgc)
            integer(ip), intent(in)  :: lshsgc
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine shsgci
    end interface

contains

    pure function get_lshses(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)              :: return_value

        ! Local variables
        integer(ip)         :: l1, l2
        type(SpherepackUtility) :: sphere_aux


        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = (l1*l2*(2*nlat-l1+1))/2+nlon+15

    end function get_lshses

    pure function get_lwork_shses(nlat, nlon) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables
        integer(ip)         :: l1, l2
        type(SpherepackUtility) :: sphere_aux

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = 5*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2

    end function get_lwork_shses

    pure function get_ldwork_shses(nlat) &
        result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value

        return_value = nlat + 1

    end function get_ldwork_shses

    pure function get_legendre_workspace_size_shses(nlat, nlon, nt, ityp) result (return_value)

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

    end function get_legendre_workspace_size_shses


    pure function get_lshsgs(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackUtility) :: sphere_aux


        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15


    end function get_lshsgs


    pure function get_lwork_shsgs(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in) :: nlat
        integer(ip)              :: return_value


        return_value = 4*nlat*(nlat+2)+2

    end function get_lwork_shsgs


    pure function get_ldwork_shsgs(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value


        return_value = nlat*(nlat+4)

    end function get_ldwork_shsgs


    pure function get_legendre_workspace_size_shsgs(nlat, nlon, nt, isym) result (return_value)

        ! Dummy arguments
        integer(ip),           intent(in) :: nlat
        integer(ip),           intent(in) :: nlon
        integer(ip), optional, intent(in) :: nt
        integer(ip), optional, intent(in) :: isym
        integer(ip)                        :: return_value

        ! Local variables
        integer(ip) :: nt_op, isym_op, l2

        !  Address optional arguments
        if (present(nt)) then
            nt_op = nt
        else
            nt_op = 1
        end if

        if (present(isym)) then
            isym_op = isym
        else
            isym_op = 0
        end if

        !
        !  Compute workspace size
        !
        select case (isym)
            case (0)
                ! Set workspace size
                return_value = nlat*nlon*(nt_op+1)
            case default
                ! Compute parity
                select case (mod(nlat, 2))
                    case (0)
                        l2 = nlat/2
                    case default
                        l2 = (nlat + 1)/2
                end select
                ! Set workspace size
                return_value = l2*nlon*(nt_op+1)
        end select

    end function get_legendre_workspace_size_shsgs


end module scalar_synthesis_routines
