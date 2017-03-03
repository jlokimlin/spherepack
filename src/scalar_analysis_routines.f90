module scalar_analysis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_SpherepackAux, only: &
        SpherepackAux

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: shaec, shaeci
    public :: shaes, shaesi
    public :: shagc, shagci
    public :: shags, shagsi
    public :: ShaesAux, ShagsAux

    ! Parameters confined to the module
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: HALF = 0.5_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: TWO = 2.0_wp
    real(wp), parameter :: FOUR = 4.0_wp

    type, public :: ShaesAux
    contains
        ! Type-bound procedures
        procedure, nopass :: shaes
        procedure, nopass :: shaesi
        procedure, nopass :: get_lshaes
        procedure, nopass :: get_lwork &
            => get_lwork_shaes
        procedure, nopass :: get_ldwork &
            => get_ldwork_shaes
        procedure, nopass :: get_legendre_workspace_size &
            => get_legendre_workspace_size_shaes
    end type ShaesAux
    
    type, public :: ShagsAux
    contains
        ! Type-bound procedures
        procedure, nopass :: shags
        procedure, nopass :: shagsi
        procedure, nopass :: get_lshags
        procedure, nopass :: get_lwork &
            => get_lwork_shags
        procedure, nopass :: get_ldwork &
            => get_ldwork_shags
        procedure, nopass :: get_legendre_workspace_size &
            => get_legendre_workspace_size_shags
    end type ShagsAux

    ! Declare interfaces for submodule implementation
    interface
        module subroutine shaec(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshaec, lshaec, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshaec(lshaec)
            integer(ip), intent(in)   :: lshaec
            real(wp),    intent(out)  :: work(lwork)
            integer(ip), intent(in)   :: lwork
            integer(ip), intent(out)  :: ierror
        end subroutine shaec

        module subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, &
            mdab, ndab, wshaes, lshaes, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshaes(lshaes)
            integer(ip), intent(in)   :: lshaes
            real(wp),    intent(out)  :: work(lwork)
            integer(ip), intent(in)   :: lwork
            integer(ip), intent(out)  :: ierror
        end subroutine shaes

        module subroutine shagc(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshagc, lshagc, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: isym
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshagc(lshagc)
            integer(ip), intent(in)   :: lshagc
            real(wp),    intent(out)  :: work(lwork)
            integer(ip), intent(in)   :: lwork
            integer(ip), intent(out)  :: ierror
        end subroutine shagc

        module subroutine shags(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshags, lshags, work, lwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            integer(ip), intent(in)   :: mode
            integer(ip), intent(in)   :: nt
            real(wp),    intent(in)   :: g(idg, jdg, nt)
            integer(ip), intent(in)   :: idg
            integer(ip), intent(in)   :: jdg
            real(wp),    intent(out)  :: a(mdab, ndab, nt)
            real(wp),    intent(out)  :: b(mdab, ndab, nt)
            integer(ip), intent(in)   :: mdab
            integer(ip), intent(in)   :: ndab
            real(wp),    intent(in)   :: wshags(lshags)
            integer(ip), intent(in)   :: lshags
            real(wp),    intent(out)  :: work(lwork)
            integer(ip), intent(in)   :: lwork
            integer(ip), intent(out)  :: ierror
        end subroutine shags

        module subroutine shaeci(nlat, nlon, wshaec, lshaec, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshaec(lshaec)
            integer(ip), intent(in)  :: lshaec
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine shaeci

        module subroutine shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, &
            ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)     :: nlat
            integer(ip), intent(in)     :: nlon
            real(wp),    intent(out)    :: wshaes(lshaes)
            integer(ip), intent(in)     :: lshaes
            real(wp),    intent(out)    :: work(lwork)
            integer(ip), intent(in)     :: lwork
            real(wp),    intent(out)    :: dwork(ldwork)
            integer(ip), intent(in)     :: ldwork
            integer(ip), intent(out)    :: ierror
        end subroutine shaesi

        module subroutine shagci(nlat, nlon, wshagc, lshagc, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)  :: nlat
            integer(ip), intent(in)  :: nlon
            real(wp),    intent(out) :: wshagc(lshagc)
            integer(ip), intent(in)  :: lshagc
            real(wp),    intent(out) :: dwork(ldwork)
            integer(ip), intent(in)  :: ldwork
            integer(ip), intent(out) :: ierror
        end subroutine shagci

        module subroutine shagsi(nlat, nlon, wshags, lshags, work, lwork, dwork, ldwork, ierror)

            ! Dummy arguments
            integer(ip), intent(in)   :: nlat
            integer(ip), intent(in)   :: nlon
            real(wp),    intent(out)  :: wshags(lshags)
            integer(ip), intent(in)   :: lshags
            real(wp),    intent(out)  :: work(lwork)
            integer(ip), intent(in)   :: lwork
            real(wp),    intent(out)  :: dwork(ldwork)
            integer(ip), intent(in)   :: ldwork
            integer(ip), intent(out)  :: ierror
        end subroutine shagsi
    end interface

contains

    pure function get_lshaes(nlat, nlon) result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables
        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = ( l1 * l2 * (2*nlat-l1+1) )/2 + nlon+15

    end function get_lshaes

    pure function get_lwork_shaes(nlat, nlon) result (return_value)

        ! Dummy arguments
        integer(ip), intent(in) :: nlat
        integer(ip), intent(in) :: nlon
        integer(ip)             :: return_value

        ! Local variables
        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = 5*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2

    end function get_lwork_shaes


    pure function get_ldwork_shaes(nlat) result (return_value)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip)              :: return_value

        return_value = nlat + 1

    end function get_ldwork_shaes


    pure function get_legendre_workspace_size_shaes(nlat, nlon, nt, ityp) result (return_value)

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

    end function get_legendre_workspace_size_shaes


    pure function get_lshags(nlat, nlon) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip)               :: return_value

        ! Local variables

        integer(ip)         :: l1, l2
        type(SpherepackAux) :: sphere_aux


        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

    end function get_lshags


    pure function get_lwork_shags(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in) :: nlat
        integer(ip)              :: return_value


        return_value = 4*nlat*(nlat+2)+2

    end function get_lwork_shags


    pure function get_ldwork_shags(nlat) result (return_value)

        ! Dummy arguments

        integer(ip), intent(in)  :: nlat
        integer(ip)               :: return_value


        return_value = nlat*(nlat+4)

    end function get_ldwork_shags


    pure function get_legendre_workspace_size_shags(nlat, nlon, nt, isym) result (return_value)

        ! Dummy arguments

        integer(ip),           intent(in) :: nlat
        integer(ip),           intent(in) :: nlon
        integer(ip), optional, intent(in) :: nt
        integer(ip), optional, intent(in) :: isym
        integer(ip)                        :: return_value

        ! Local variables

        integer(ip) :: nt_op, isym_op, l2


        !
        !  Address optional arguments
        !
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

    end function get_legendre_workspace_size_shags


end module scalar_analysis_routines
