module scalar_analysis_routines

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    use module_gaqd, only: &
        gaqd

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    !public :: shaec
    !public :: shaeci
    public :: shaes
    public :: shaesi
    public :: ShaesAux
    !public :: shagc
    !public :: shagci
    public :: shags
    public :: shagsi
    public :: ShagsAux

    ! Declare derived data type
    type, public :: ShaesAux
        !-----------------------------------------
        ! Type components
        !-----------------------------------------
    contains
        !-----------------------------------------
        ! Type-bound procedures
        !-----------------------------------------
        procedure, nopass :: shaes
        procedure, nopass :: shaesi
        procedure, nopass :: get_lshaes
        procedure, nopass :: get_lwork &
            => get_lwork_shaes
        procedure, nopass :: get_ldwork &
            => get_ldwork_shaes
        procedure, nopass :: get_legendre_workspace_size &
            => get_legendre_workspace_size_shaes
        !-----------------------------------------
    end type ShaesAux

    ! Declare derived data type
    type, public :: ShagsAux
        !-----------------------------------------
        ! Type components
        !-----------------------------------------
    contains
        !-----------------------------------------
        ! Type-bound procedures
        !-----------------------------------------
        procedure, nopass :: shags
        procedure, nopass :: shagsi
        procedure, nopass :: get_lshags
        procedure, nopass :: get_lwork &
            => get_lwork_shags
        procedure, nopass :: get_ldwork &
            => get_ldwork_shags
        procedure, nopass :: get_legendre_workspace_size &
            => get_legendre_workspace_size_shags
        !-----------------------------------------
    end type ShagsAux

    interface

        module subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, &
            mdab, ndab, wshaes, lshaes, work, lwork, ierror)
            !----------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            integer (ip), intent (in)     :: isym
            integer (ip), intent (in)     :: nt
            real (wp),    intent (in)     :: g(idg, jdg, nt)
            integer (ip), intent (in)     :: idg
            integer (ip), intent (in)     :: jdg
            real (wp),    intent (out)    :: a(mdab, ndab, nt)
            real (wp),    intent (out)    :: b(mdab, ndab, nt)
            integer (ip), intent (in)     :: mdab
            integer (ip), intent (in)     :: ndab
            real (wp),    intent (in out) :: wshaes(lshaes)
            integer (ip), intent (in)     :: lshaes
            real (wp),    intent (in out) :: work(lwork)
            integer (ip), intent (in)     :: lwork
            integer (ip), intent (out)    :: ierror
            !----------------------------------------------------------
        end subroutine shaes

        module subroutine shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, &
            ldwork, ierror)
            !
            ! Remarks:
            !
            ! size(wshaes) = (l*(l+1)*imid)/2+nlon+15
            ! size(work) = 5*l*imid + 3*((l-3)*l+2)/2
            !
            !----------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            real (wp),    intent (out)    :: wshaes(lshaes)
            integer (ip), intent (in)     :: lshaes
            real (wp),    intent (out)    :: work(lwork)
            integer (ip), intent (in)     :: lwork
            real (wp),    intent (out)    :: dwork(ldwork)
            integer (ip), intent (in)     :: ldwork
            integer (ip), intent (out)    :: ierror
            !----------------------------------------------------------
        end subroutine shaesi

        module subroutine shags(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, &
            wshags, lshags, work, lwork, ierror)
            !
            ! Purpose:
            !
            ! Performs the spherical harmonic analysis on
            ! a gaussian grid on the array(s) in g and returns the coefficients
            ! in array(s) a, b. the necessary legendre polynomials are fully
            ! stored in this version.
            !
            !----------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            integer (ip), intent (in)     :: mode
            integer (ip), intent (in)     :: nt
            real (wp),    intent (in)     :: g(idg, jdg, nt)
            integer (ip), intent (in)     :: idg
            integer (ip), intent (in)     :: jdg
            real (wp),    intent (out)    :: a(mdab, ndab, nt)
            real (wp),    intent (out)    :: b(mdab, ndab, nt)
            integer (ip), intent (in)     :: mdab
            integer (ip), intent (in)     :: ndab
            real (wp),    intent (in out) :: wshags(lshags)
            integer (ip), intent (in)     :: lshags
            real (wp),    intent (in out) :: work(lwork)
            integer (ip), intent (in)     :: lwork
            integer (ip), intent (out)    :: ierror
            !----------------------------------------------------------
        end subroutine shags

        module subroutine shagsi(nlat, nlon, wshags, lshags, work, lwork, dwork, ldwork, ierror)
            !
            !     Remark:
            !
            !     this subroutine must be called before calling shags or shsgs with
            !     fixed nlat, nlon. it precomputes the gaussian weights, points
            !     and all necessary legendre polys and stores them in wshags.
            !     these quantities must be preserved when calling shags or shsgs
            !     repeatedly with fixed nlat, nlon.  dwork must be of length at
            !     least nlat*(nlat+4) in the routine calling shagsi.  This is
            !     not checked.  undetectable errors will result if dwork is
            !     smaller than nlat*(nlat+4).
            !
            !----------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            real (wp),    intent (in out) :: wshags(lshags)
            integer (ip), intent (in)     :: lshags
            real (wp),    intent (in out) :: work(lwork)
            integer (ip), intent (in)     :: lwork
            real (wp),    intent (in out) :: dwork(ldwork)
            integer (ip), intent (in)     :: ldwork
            integer (ip), intent (out)    :: ierror
            !----------------------------------------------------------
        end subroutine shagsi

    end interface

contains


    pure function get_lshaes(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)         :: l1, l2
        type (SpherepackAux) :: sphere_aux
        !----------------------------------------------------------------------

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = ( l1 * l2 * (2*nlat-l1+1) )/2 + nlon+15

    end function get_lshaes


    pure function get_lwork_shaes(nlat, nlon) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)         :: l1, l2
        type (SpherepackAux) :: sphere_aux
        !----------------------------------------------------------------------

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = 5*nlat*l2+3*((l1-2)*(2*nlat-l1-1))/2

    end function get_lwork_shaes


    pure function get_ldwork_shaes(nlat) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip)               :: return_value
        !----------------------------------------------------------------------

        return_value = nlat + 1

    end function get_ldwork_shaes


    pure function get_legendre_workspace_size_shaes(nlat, nlon, nt, ityp) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip),           intent (in) :: nlat
        integer (ip),           intent (in) :: nlon
        integer (ip), optional, intent (in) :: nt
        integer (ip), optional, intent (in) :: ityp
        integer (ip)                        :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: nt_op, ityp_op, l2
        !----------------------------------------------------------------------

        !
        !==> Address optional arguments
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
        !==> Compute workspace size
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
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip)               :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)         :: l1, l2
        type (SpherepackAux) :: sphere_aux
        !----------------------------------------------------------------------

        call sphere_aux%compute_parity(nlat, nlon, l1, l2)

        return_value = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15

    end function get_lshags


    pure function get_lwork_shags(nlat) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat
        integer (ip)              :: return_value
        !----------------------------------------------------------------------

        return_value = 4*nlat*(nlat+2)+2

    end function get_lwork_shags


    pure function get_ldwork_shags(nlat) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip)               :: return_value
        !----------------------------------------------------------------------

        return_value = nlat*(nlat+4)

    end function get_ldwork_shags


    pure function get_legendre_workspace_size_shags(nlat, nlon, nt, isym) result (return_value)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip),           intent (in) :: nlat
        integer (ip),           intent (in) :: nlon
        integer (ip), optional, intent (in) :: nt
        integer (ip), optional, intent (in) :: isym
        integer (ip)                        :: return_value
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: nt_op, isym_op, l2
        !----------------------------------------------------------------------

        !
        !==> Address optional arguments
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
        !==> Compute workspace size
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
