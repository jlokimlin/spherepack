!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK version 3.2                   *
!     *                                                               *
!     *       A Package of Fortran Subroutines and Programs           *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! ... file type_SpherepackAux.f90
!
!     This file must be loaded with all main program files
!     in spherepack. It includes undocumented subroutines
!     called by some or all of main programs
!
module type_SpherepackAux

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: SpherepackAux

    type, public :: SpherepackAux
        !----------------------------------------------------------------------
        ! Type components
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Type-bound procedures
        !----------------------------------------------------------------------
        procedure, nopass :: alin
        procedure, nopass :: alinit
        procedure, nopass :: compute_parity
        procedure, nopass :: dnlfk
        procedure, nopass :: dnlft
        procedure, nopass :: dnlftd
        procedure, nopass :: dnzfk
        procedure, nopass :: dnzft
        procedure, nopass :: dvbk
        procedure, nopass :: dvbt
        procedure, nopass :: dvtk
        procedure, nopass :: dvtt
        procedure, nopass :: dwbk
        procedure, nopass :: dwbt
        procedure, nopass :: dwtk
        procedure, nopass :: dwtt
        procedure, nopass :: dzvk
        procedure, nopass :: dzvt
        procedure, nopass :: dzwk
        procedure, nopass :: dzwt
        procedure, nopass :: legin
        procedure, nopass :: rabcp
        procedure, nopass :: rabcv
        procedure, nopass :: rabcw
        procedure, nopass :: sea1
        procedure, nopass :: ses1
        procedure, nopass :: vbgint
        procedure, nopass :: vbin
        procedure, nopass :: vbinit
        procedure, nopass :: vtgint
        procedure, nopass :: vtinit
        procedure, nopass :: wbgint
        procedure, nopass :: wbin
        procedure, nopass :: wbinit
        procedure, nopass :: wtgint
        procedure, nopass :: wtinit
        procedure, nopass :: zfin
        procedure, nopass :: zfinit
        procedure, nopass :: zvin
        procedure, nopass :: zvinit
        procedure, nopass :: zwin
        procedure, nopass :: zwinit
        !----------------------------------------------------------------------
    end type SpherepackAux


contains


    pure subroutine compute_parity(nlat, nlon, l1, l2)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip), intent (out) :: l1
        integer (ip), intent (out) :: l2
        !----------------------------------------------------------------------

        ! Compute parity in nlon
        select case (mod(nlon, 2))
            case (0)
                l1 = min(nlat, (nlon+2)/2)
            case default
                l1 = min(nlat, (nlon+1)/2)
        end select

        ! Compute parity in nlat
        select case (mod(nlat, 2))
            case (0)
                l2 = nlat/2
            case default
                l2 = (nlat + 1)/2
        end select

    end subroutine compute_parity


    pure subroutine dnlfk(m, n, cp)
        !
        !     cp requires n/2+1 real locations
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: cp(n/2+1)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)         :: i, l, ma, nex,  nmms2
        real (wp)            :: a1, b1, c1, t1, t2
        real (wp)            :: fk, cp2, pm1
        real (wp)            :: fden, fnmh, fnum, fnnp1, fnmsq
        real (wp), parameter :: SC10=1024.0_wp
        real (wp), parameter :: SC20=SC10**2
        real (wp), parameter :: SC40=SC20**2
        !----------------------------------------------------------------------

        ma = abs(m)

        if (ma > n) then
            cp(1) = 0.0_wp
            return
        end if

        if (n < 1) then
            cp(1) = sqrt(2.0_wp)
        else if (n == 1) then
            if (ma /= 0) then
                cp(1) = sqrt(0.75_wp)
                if (m == -1) then
                    cp(1) = -cp(1)
                end if
            else
                cp(1) = sqrt(1.5_wp)
            end if
        else
            if (mod(n+ma, 2) /= 0) then
                nmms2 = (n-ma-1)/2
                fnum = real(n+ma+2, kind=wp)
                fnmh = real(n-ma+2, kind=wp)
                pm1 = -1.0_wp
            else
                nmms2 = (n-ma)/2
                fnum = real(n+ma+1, kind=wp)
                fnmh = real(n-ma+1, kind=wp)
                pm1 = 1.0_wp
            end if

            t1 = 1.0_wp/SC20
            nex = 20
            fden = 2.0_wp

            if (nmms2 >= 1) then
                do i=1, nmms2
                    t1 = fnum*t1/fden
                    if (t1 > SC20) then
                        t1 = t1/SC40
                        nex = nex+40
                    end if
                    fnum = fnum+2.0_wp
                    fden = fden+2.0_wp
                end do
            end if

            t1 = t1/2.0**(n-1-nex)

            if (mod(ma/2, 2) /= 0) t1 = -t1

            t2 = 1.0_wp

            if (ma /= 0) then
                do i=1, ma
                    t2 = fnmh*t2/(fnmh+pm1)
                    fnmh = fnmh+2.0_wp
                end do
            end if

            cp2 = t1*sqrt((real(n, kind=wp)+0.5_wp)*t2)
            fnnp1 = real(n*(n+1), kind=wp)
            fnmsq = fnnp1 - 2.0_wp * real(ma**2, kind=wp)
            l = (n+1)/2

            if (mod(n, 2) == 0 .and. mod(ma, 2) == 0) l = l+1

            cp(l) = cp2

            if ((m < 0) .and. (mod(ma, 2) /= 0)) cp(l) = -cp(l)

            if (l <= 1) return

            fk = real(n, kind=wp)
            a1 = (fk-2.0_wp)*(fk-1.0_wp)-fnnp1
            b1 = 2.0_wp*(fk**2-fnmsq)
            cp(l-1) = b1*cp(l)/a1

            l = l - 1

            do while (l > 1)
                fk = fk-2.0_wp
                a1 = (fk-2.0_wp)*(fk-1.0_wp)-fnnp1
                b1 = -2.0_wp*(fk**2-fnmsq)
                c1 = (fk+1.0_wp)*(fk+2.0_wp)-fnnp1
                cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
                l=l-1
            end do
        end if

    end subroutine dnlfk


    subroutine dnlftd (m, n, theta, cp, pb)
        !
        ! Purpose:
        !
        ! Computes the derivative of pmn(theta) with respect to theta
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (in)  :: theta
        real (wp),    intent (out) :: cp(*)
        real (wp),    intent (out) :: pb
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) ::  k, kdo
        real (wp)    :: cos2t, sin2t, cost, sint, temp
        !----------------------------------------------------------------------

        cos2t = cos(2.0_wp*theta)
        sin2t = sin(2.0_wp*theta)

        if (mod(n, 2) <= 0) then
            if (mod(abs(m), 2) <= 0) then
                !
                !==> n even, m even
                !
                kdo=n/2
                pb = 0.0_wp

                if (n == 0) return

                cost = cos2t
                sint = sin2t
                do k=1, kdo
                    !     pb = pb+cp(k+1)*cos(2*k*theta)
                    pb = pb-2.0_wp*real(k, kind=wp)*cp(k+1)*sint
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            else
                !
                !==> n even, m odd
                !
                kdo = n/2
                pb = 0.0_wp
                cost = cos2t
                sint = sin2t
                do k=1, kdo
                    !     pb = pb+cp(k)*sin(2*k*theta)
                    pb = pb+2.0_wp*real(k, kind=wp)*cp(k)*cost
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            end if
        else
            if (mod(abs(m), 2) <= 0) then
                !
                !==> n odd, m even
                !
                kdo = (n+1)/2
                pb = 0.0_wp
                cost = cos(theta)
                sint = sin(theta)
                do k=1, kdo
                    !     pb = pb+cp(k)*cos((2*k-1)*theta)
                    pb = pb-(2.0_wp*real(k-1, kind=wp))*cp(k)*sint
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            else
                !
                !==> n odd, m odd
                !
                kdo = (n+1)/2
                pb = 0.0_wp
                cost = cos(theta)
                sint = sin(theta)
                do k=1, kdo
                    !     pb = pb+cp(k)*sin((2*k-1)*theta)
                    pb = pb+(2.0_wp*real(k-1, kind=wp))*cp(k)*cost
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            end if
        end if

    end subroutine dnlftd


    pure subroutine dnlft(m, n, theta, cp, pb)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (in)  :: theta
        real (wp),    intent (out) :: cp(*)
        real (wp),    intent (out) :: pb
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) ::  k, kdo
        real (wp)    :: temp, cos2t, cost, sin2t, sint
        !----------------------------------------------------------------------

        cos2t = cos(2.0_wp * theta)
        sin2t = sin(2.0_wp * theta)

        if (mod(n, 2) <= 0) then
            if (mod(m, 2) <= 0) then
                !
                !==>  n even, m even
                !
                kdo = n/2
                pb = 0.5_wp * cp(1)

                if (n == 0) return

                cost = cos2t
                sint = sin2t

                do k=1, kdo
                    pb = pb+cp(k+1)*cost
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            else
                !
                !==> n even, m odd
                !
                kdo = n/2
                pb = 0.0_wp
                cost = cos2t
                sint = sin2t

                do k=1, kdo
                    pb = pb+cp(k)*sint
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            end if
        else
            if (mod(m, 2) <= 0) then
                !
                !==> n odd, m even
                !
                kdo = (n+1)/2
                pb = 0.0_wp
                cost = cos(theta)
                sint = sin(theta)

                do k=1, kdo
                    pb = pb+cp(k)*cost
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            else
                !
                !==>  n odd, m odd
                !
                kdo = (n+1)/2
                pb = 0.0_wp
                cost = cos(theta)
                sint = sin(theta)

                do k=1, kdo
                    pb = pb+cp(k)*sint
                    temp = cos2t*cost-sin2t*sint
                    sint = sin2t*cost+cos2t*sint
                    cost = temp
                end do
            end if
        end if

    end subroutine dnlft


    subroutine legin(mode, l, nlat, m, w, pmn, km)
        !
        ! Purpose:
        !
        ! Computes legendre polynomials for n=m, ..., l-1
        ! and  i=1, ..., late (late=((nlat+mod(nlat, 2))/2)gaussian grid
        ! in pmn(n+1, i, km) using swarztrauber's recursion formula.
        ! the vector w contains quantities precomputed in shigc.
        ! legin must be called in the order m=0, 1, ..., l-1
        ! (e.g., if m=10 is sought it must be preceded by calls with
        ! m=0, 1, 2, ..., 9 in that order)
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: mode
        integer (ip), intent (in)  :: l
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: m
        real (wp),    intent (out) :: w(*)
        real (wp),    intent (out) :: pmn(*)
        integer (ip), intent (out) :: km
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: late
        integer (ip) :: workspace_indices(5)
        !----------------------------------------------------------------------

        !
        !==> set size of pole to equator gaussian grid
        !
        late = (nlat+mod(nlat,2))/2

        !
        !==> partition w (set pointers for p0n, p1n, abel, bbel, cbel, pmn)
        !
        workspace_indices = get_workspace_indices(l, late, nlat)

        associate( i => workspace_indices )

            call legin1(mode, l, nlat, late, m, &
                w(i(1)), w(i(2)), w(i(3)), w(i(4)), w(i(5)), pmn, km)

        end associate


    contains


        subroutine legin1(mode, l, nlat, late, m, p0n, p1n, abel, bbel, cbel, &
            pmn, km)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)  :: mode
            integer (ip), intent (in)  :: l
            integer (ip), intent (in)  :: nlat
            integer (ip), intent (in)  :: late
            integer (ip), intent (in)  :: m
            real (wp),    intent (out) :: p0n(nlat, late)
            real (wp),    intent (out) :: p1n(nlat, late)
            real (wp),    intent (out) :: abel(*)
            real (wp),    intent (out) :: bbel(*)
            real (wp),    intent (out) :: cbel(*)
            real (wp),    intent (out) :: pmn(nlat, late, 3)
            integer (ip), intent (out) :: km
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip)       :: i, n, ms, np1, imn, kmt, ninc
            integer (ip), save :: column_indices(0:2) = [1, 2, 3]
            !----------------------------------------------------------------------

            !     set do loop indices for full or half sphere
            ms = m+1
            ninc = 1
            select case (mode)
                case (1)
                    !     only compute pmn for n-m odd
                    ms = m+2
                    ninc = 2
                case (2)
                    !     only compute pmn for n-m even
                    ms = m+1
                    ninc = 2
            end select

            associate( &
                km0 => column_indices(0), &
                km1 => column_indices(1), &
                km2 => column_indices(2) &
                )

                if (m > 1) then
                    do np1=ms, nlat, ninc
                        n = np1-1
                        imn = indx(m, n)
                        if (n >= l) imn = imndx(l, m, n)
                        do i=1, late
                            pmn(np1, i, km0) = abel(imn)*pmn(n-1, i, km2) &
                                +bbel(imn)*pmn(n-1, i, km0) &
                                -cbel(imn)*pmn(np1, i, km2)
                        end do
                    end do
                else if (m == 0) then
                    do np1=ms, nlat, ninc
                        do i=1, late
                            pmn(np1, i, km0) = p0n(np1, i)
                        end do
                    end do
                else if (m == 1) then
                    do np1=ms, nlat, ninc
                        do i=1, late
                            pmn(np1, i, km0) = p1n(np1, i)
                        end do
                    end do
                end if

                !
                !==>  Permute column indices
                !     km0, km1, km2 store m, m-1, m-2 columns
                kmt = km0
                km0 = km2
                km2 = km1
                km1 = kmt
            end associate

            !     set current m index in output param km
            km = kmt

        end subroutine legin1


        pure function get_workspace_indices(l, late, nlat) &
            result (return_value)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in) :: l
            integer (ip), intent (in) :: late
            integer (ip), intent (in) :: nlat
            integer (ip)              :: return_value(5)
            !----------------------------------------------------------------------

            associate( i => return_value )

                i(1) = 1+nlat
                i(2) = i(1)+nlat*late
                i(3) = i(2)+nlat*late
                i(4) = i(3)+(2*nlat-l)*(l-1)/2
                i(5) = i(4)+(2*nlat-l)*(l-1)/2

            end associate

        end function get_workspace_indices


        pure function indx(m, n) result (return_value)
            !
            ! Purpose:
            !
            !     index function used in storing triangular
            !     arrays for recursion coefficients (functions of (m, n))
            !     for 2 <= m <= n-1 and 2 <= n <= l-1
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in) :: m
            integer (ip), intent (in) :: n
            integer (ip)              :: return_value
            !----------------------------------------------------------------------


            return_value = (n-1)*(n-2)/2+m-1


        end function indx

        pure function imndx(l, m, n) result (return_value)
            !
            ! Purpose:
            !
            !     index function used in storing triangular
            !     arrays for recursion coefficients (functions of (m, n))
            !     for l <= n <= nlat and 2 <= m <= l
            !
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in) :: l
            integer (ip), intent (in) :: m
            integer (ip), intent (in) :: n
            integer (ip)              :: return_value
            !----------------------------------------------------------------------

            return_value = l*(l-1)/2+(n-l-1)*(l-1)+m-1

        end function imndx

    end subroutine legin



    subroutine zfin(isym, nlat, nlon, m, z, i3, wzfin)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: isym
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        integer (ip), intent (in)     :: m
        real (wp),    intent (in out) :: z(*)
        integer (ip), intent (in out) :: i3
        real (wp),    intent (in out) :: wzfin(*)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: imid
        integer (ip) :: workspace(4)
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        !
        !==> The length of wzfin is 2*lim+3*labc
        !
        workspace = get_workspace_indices(nlat, nlon, imid)

        associate( &
            iw1 => workspace(1), &
            iw2 => workspace(2), &
            iw3 => workspace(3), &
            iw4 => workspace(4) &
            )

            call zfin1(isym, nlat, m, z, imid, i3, wzfin, &
                wzfin(iw1), wzfin(iw2), wzfin(iw3), wzfin(iw4))

        end associate


    contains

        subroutine zfin1(isym, nlat, m, z, imid, i3, zz, z1, a, b, c)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: isym
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: m
            real (wp),    intent (in out) :: z(imid, nlat, 3)
            integer (ip), intent (in)     :: imid
            integer (ip), intent (in out) :: i3
            real (wp),    intent (in out) :: zz(imid,*)
            real (wp),    intent (in out) :: z1(imid,*)
            real (wp),    intent (in out) :: a(*)
            real (wp),    intent (in out) :: b(*)
            real (wp),    intent (in out) :: c(*)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), save :: i1, i2
            integer (ip)       :: ns, np1, nstp, itemp, nstrt
            !----------------------------------------------------------------------

            itemp = i1
            i1 = i2
            i2 = i3
            i3 = itemp

            if (m < 1) then
                i1 = 1
                i2 = 2
                i3 = 3
                z(:,1:nlat, i3) = zz(:,1:nlat)
            else if (m == 1) then
                z(:,2:nlat, i3) = z1(:,2:nlat)
            else
                ns = ((m-2)*(nlat+nlat-m-1))/2+1

                if (isym /= 1) then
                    z(:,m+1,i3) = a(ns)*z(:,m-1,i1)-c(ns)*z(:,m+1,i1)
                end if

                if (m == nlat-1) return

                if (isym /= 2) then
                    ns = ns+1
                    z(:, m+2, i3) = a(ns)*z(:, m, i1) -c(ns)*z(:, m+2, i1)
                end if

                nstrt = m+3

                if (isym == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2

                if (isym == 0) nstp = 1

                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    z(:,np1,i3) = &
                        a(ns)*z(:,np1-2,i1)+b(ns)*z(:,np1-2,i3)-c(ns)*z(:,np1,i1)
                end do
            end if

        end subroutine zfin1



        pure function get_workspace_indices(nlat, nlon, imid) &
            result (return_value)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip), intent (in) :: imid
            integer (ip)              :: return_value(4)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip) :: mmax, lim, labc
            !----------------------------------------------------------------------

            associate( i => return_value )

                mmax = min(nlat, nlon/2+1)
                lim = nlat*imid
                labc = ((mmax-2)*(2*nlat-mmax-1))/2
                i(1) = lim+1
                i(2) = i(1)+lim
                i(3) = i(2)+labc
                i(4) = i(3)+labc

            end associate

        end function get_workspace_indices

    end subroutine zfin


    subroutine zfinit(nlat, nlon, wzfin, dwork)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        real (wp),    intent (in out) :: wzfin(*)
        real (wp),    intent (in out) :: dwork(nlat+2)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid
        integer (ip) :: iw1, iw2
        !----------------------------------------------------------------------

        imid = (nlat+1)/2

        !
        ! Remarks:
        !
        !     the length of wzfin is 3*((l-3)*l+2)/2 + 2*l*imid
        !     the length of dwork is nlat+2
        !

        ! Set workspace indices
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+1

        call zfini1(nlat, nlon, imid, wzfin, wzfin(iw1), dwork, dwork(iw2))

    contains

        subroutine zfini1(nlat, nlon, imid, z, abc, cz, work)
            !
            !     Remarks:
            !
            !     abc must have 3*((mmax-2)*(2*nlat-mmax-1))/2 locations
            !     where mmax = min(nlat, nlon/2+1)
            !     cz and work must each have nlat+1 locations
            !
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            integer (ip), intent (in)     :: imid
            real (wp),    intent (in out) :: z(imid, nlat, 2)
            real (wp),    intent (in out) :: abc(*)
            real (wp),    intent (in out) :: cz(nlat+1)
            real (wp),    intent (in out) :: work(nlat+1)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip)         :: i, m, n, mp1, np1
            real (wp)            :: dt, th, zh
            !----------------------------------------------------------------------

            dt = PI/(nlat-1)

            do mp1=1, 2
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dnzfk(nlat, m, n, cz, work)
                    do i=1, imid
                        th = real(i-1, kind=wp)*dt
                        call dnzft(nlat, m, n, th, cz, zh)
                        z(i, np1, mp1) = zh
                    end do
                    z(1, np1, mp1) = 0.5_wp*z(1, np1, mp1)
                end do
            end do

            call rabcp(nlat, nlon, abc)

        end subroutine zfini1

    end subroutine zfinit



    subroutine dnzfk(nlat, m, n, cz, work)
        !
        ! Purpose:
        !
        !     Computes the coefficients in the trigonometric
        !     expansion of the z functions that are used in spherical
        !     harmonic analysis.
        !
        !
        ! Remark:
        !
        !    cz and work must both have nlat/2+1 locations
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: n
        real (wp),    intent (in out) :: cz(nlat/2 + 1)
        real (wp),    intent (in out) :: work(nlat/2 + 1)
        !----------------------------------------------------------------------
        integer (ip) :: i, k, lc, kp1, kdo, idx
        real (wp)    :: summation, sc1, t1, t2
        !----------------------------------------------------------------------

        lc = (nlat+1)/2
        sc1 = 2.0_wp/(nlat-1)

        call dnlfk(m, n, work)

        if (mod(n,2) <= 0) then
            if (mod(m,2) <= 0) then
                kdo = n/2+1
                do idx=1, lc
                    i = 2*idx-2
                    summation = work(1)/(1.0_wp-real(i**2, kind=wp))
                    if (kdo >= 2) then
                        do kp1=2, kdo
                            k = kp1-1
                            t1 = 1.0_wp-real((2*k+i)**2, kind=wp)
                            t2 = 1.0_wp-real((2*k-i)**2, kind=wp)
                            summation = summation+work(kp1)*(t1+t2)/(t1*t2)
                        end do
                    end if
                    cz(idx) = sc1*summation
                end do
            else
                kdo = n/2
                do idx=1, lc
                    i = 2*idx-2
                    summation = 0.0_wp
                    do k=1, kdo
                        t1 = 1.0_wp-real((2*k+i)**2, kind=wp)
                        t2 = 1.0_wp-real((2*k-i)**2, kind=wp)
                        summation=summation+work(k)*(t1-t2)/(t1*t2)
                    end do
                    cz(idx) = sc1*summation
                end do
            end if
        else
            if (mod(m,2) <= 0) then
                !
                !==>  n odd, m even
                !
                kdo = (n+1)/2
                do idx=1, lc
                    i = 2*idx-1
                    summation = 0.0_wp
                    do k=1, kdo
                        t1 = 1.0_wp-real((2*k-1+i)**2, kind=wp)
                        t2 = 1.0_wp-real((2*k-1-i)**2, kind=wp)
                        summation=summation+work(k)*(t1+t2)/(t1*t2)
                    end do
                    cz(idx)=sc1*summation
                end do
            else
                kdo = (n+1)/2
                do idx=1, lc
                    i = 2*idx-3
                    summation=0.0_wp
                    do k=1, kdo
                        t1 = 1.0_wp-real((2*k-1+i)**2, kind=wp)
                        t2 = 1.0_wp-real((2*k-1-i)**2, kind=wp)
                        summation=summation+work(k)*(t1-t2)/(t1*t2)
                    end do
                    cz(idx)=sc1*summation
                end do
            end if
        end if

    end subroutine dnzfk



    subroutine dnzft(nlat, m, n, th, cz, zh)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)   :: nlat
        integer (ip), intent (in)   :: m
        integer (ip), intent (in)   :: n
        real (wp),    intent (in)   :: th
        real (wp),    intent (in)   :: cz(*)
        real (wp),    intent (out)  :: zh
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: i, k, lc, lq, ls
        real (wp)    :: cos2t, sin2t, cost, sint, temp
        !----------------------------------------------------------------------

        zh = 0.0_wp
        cos2t = cos(2.0_wp*th)
        sin2t = sin(2.0_wp*th)

        if (mod(nlat, 2) <= 0) then
            lc = nlat/2
            lq = lc-1
            if (mod(n, 2) <= 0) then
                if (mod(m, 2) <= 0) then
                    zh = 0.5_wp*cz(1)
                    cost = cos2t
                    sint = sin2t
                    do k=2, lc
                        !     zh = zh+cz(k)*cos(2*(k-1)*th)
                        zh = zh+cz(k)*cost
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                else
                    cost = cos2t
                    sint = sin2t
                    do k=1, lq
                        !     zh = zh+cz(k+1)*sin(2*k*th)
                        zh = zh+cz(k+1)*sint
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                end if
            else
                if (mod(m, 2) <= 0) then
                    zh = 0.5_wp*cz(lc)*cos((nlat-1)*th)
                    cost = cos(th)
                    sint = sin(th)
                    do k=1, lq
                        !     zh = zh+cz(k)*cos((2*k-1)*th)
                        zh = zh+cz(k)*cost
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                else
                    cost = cos(th)
                    sint = sin(th)
                    do k=1, lq
                        !     zh = zh+cz(k+1)*sin((2*k-1)*th)
                        zh = zh+cz(k+1)*sint
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                end if
            end if
        else
            lc = (nlat+1)/2
            lq = lc-1
            ls = lc-2
            if (mod(n, 2) <= 0) then
                if (mod(m, 2) <= 0) then
                    zh = 0.5_wp*(cz(1)+cz(lc)*cos(2*lq*th))
                    cost = cos2t
                    sint = sin2t
                    do k=2, lq
                        !     zh = zh+cz(k)*cos(2*(k-1)*th)
                        zh = zh+cz(k)*cost
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                else
                    cost = cos2t
                    sint = sin2t
                    do k=1, ls
                        !     zh = zh+cz(k+1)*sin(2*k*th)
                        zh = zh+cz(k+1)*sint
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                end if
            else
                if (mod(m, 2) <= 0) then
                    cost = cos(th)
                    sint = sin(th)
                    do k=1, lq
                        !     zh = zh+cz(k)*cos((2*k-1)*th)
                        zh = zh+cz(k)*cost
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                else
                    cost = cos(th)
                    sint = sin(th)
                    do k=1, lq
                        !     zh = zh+cz(k+1)*sin((2*k-1)*th)
                        zh = zh+cz(k+1)*sint
                        temp = cos2t*cost-sin2t*sint
                        sint = sin2t*cost+cos2t*sint
                        cost = temp
                    end do
                end if
            end if
        end if

    end subroutine dnzft



    subroutine alin(isym, nlat, nlon, m, p, i3, walin)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: isym
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        integer (ip), intent (in)     :: m
        real (wp),    intent (in out) :: p(*)
        integer (ip), intent (in out) :: i3
        real (wp),    intent (in out) :: walin(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid
        integer (ip) :: workspace_indices(4)
        !----------------------------------------------------------------------

        imid = (nlat+1)/2

        workspace_indices = get_workspace_indices(nlat, nlon, imid)

        associate( &
            iw1 => workspace_indices(1), &
            iw2 => workspace_indices(2), &
            iw3 => workspace_indices(3), &
            iw4 => workspace_indices(4) &
            )
            !
            !     the length of walin is ((5*l-7)*l+6)/2
            !
            call alin1(isym, nlat, m, p, imid, i3, walin, walin(iw1), walin(iw2), &
                walin(iw3), walin(iw4))

        end associate

    contains


        subroutine alin1(isym, nlat, m, p, imid, i3, pz, p1, a, b, c)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: isym
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: m
            real (wp),    intent (in out) :: p(imid,nlat,3)
            integer (ip), intent (in)     :: imid
            integer (ip), intent (in out) :: i3
            real (wp),    intent (in out) :: pz(imid,*)
            real (wp),    intent (in out) :: p1(imid,*)
            real (wp),    intent (in out) :: a(*)
            real (wp),    intent (in out) :: b(*)
            real (wp),    intent (in out) :: c(*)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip)       :: ns, np1, nstp, itemp, nstrt
            integer (ip), save :: i1, i2
            !----------------------------------------------------------------------

            itemp = i1
            i1 = i2
            i2 = i3
            i3 = itemp

            if (m < 1) then
                i1 = 1
                i2 = 2
                i3 = 3
                p(:,:, i3) = pz(:,1:nlat)
            else if (m == 1) then
                p(:,2:nlat, i3) = p1(:,2:nlat)
            else
                ns = ((m-2)*(nlat+nlat-m-1))/2+1

                if (isym /= 1) p(:, m+1, i3) = a(ns)*p(:, m-1, i1)-c(ns)*p(:, m+1, i1)

                if (m == nlat-1) return

                if (isym /= 2) then
                    ns = ns+1
                    p(:, m+2, i3) = a(ns)*p(:, m, i1)-c(ns)*p(:, m+2, i1)
                end if

                nstrt = m+3

                if (isym == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2

                if (isym == 0) nstp = 1

                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    p(:,np1,i3) = &
                        a(ns)*p(:,np1-2,i1)+b(ns)*p(:,np1-2,i3)-c(ns)*p(:,np1,i1)
                end do
            end if

        end subroutine alin1

        pure function get_workspace_indices(nlat, nlon, imid) &
            result (return_value)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip), intent (in) :: imid
            integer (ip)              :: return_value(4)
            !----------------------------------------------------------------------
            ! Local variables arguments
            !----------------------------------------------------------------------
            integer (ip) :: lim, mmax, labc
            !----------------------------------------------------------------------

            associate( i => return_value )

                lim = nlat*imid
                mmax = min(nlat, nlon/2+1)
                labc = ((mmax-2)*(2*nlat-mmax-1))/2
                i(1) = lim+1
                i(2) = i(1)+lim
                i(3) = i(2)+labc
                i(4) = i(3)+labc

            end associate

        end function get_workspace_indices

    end subroutine alin



    subroutine alinit(nlat, nlon, walin, dwork)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        real (wp),    intent (in out) :: walin(*)
        real (wp),    intent (in out) :: dwork(nlat+1)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: imid, iw1
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        !
        !     the length of walin is 3*((l-3)*l+2)/2 + 2*l*imid
        !     the length of work is nlat+1
        !
        call alini1(nlat, nlon, imid, walin, walin(iw1), dwork)

    contains

        subroutine alini1(nlat, nlon, imid, p, abc, cp)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            integer (ip), intent (in)     :: imid
            real (wp),    intent (in out) :: p(imid, nlat,2)
            real (wp),    intent (in out) :: abc(*)
            real (wp),    intent (in out) :: cp(*)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip)         :: i, m, n, mp1, np1
            real (wp)            :: dt, ph, th
            !----------------------------------------------------------------------

            dt = PI/(nlat-1)

            do mp1=1, 2
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dnlfk(m, n, cp)
                    do i=1, imid
                        th = real(i-1, kind=wp)*dt
                        call dnlft(m, n, th, cp, ph)
                        p(i, np1, mp1) = ph
                    end do
                end do
            end do

            call rabcp(nlat, nlon, abc)

        end subroutine alini1

    end subroutine alinit


    subroutine rabcp(nlat, nlon, abc)
        !
        ! Purpose:
        !
        ! Computes the coefficients in the recurrence
        ! relation for the associated legendre functions. array abc
        ! must have 3*((mmax-2)*(2*nlat-mmax-1))/2 locations.
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        real (wp),    intent (in out) :: abc(*)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: mmax, labc, iw1, iw2
        !----------------------------------------------------------------------

        ! Compute workspace indices
        mmax = min(nlat, nlon/2+1)
        labc = ((mmax-2)*(2*nlat-mmax-1))/2
        iw1 = labc+1
        iw2 = iw1+labc

        call rabcp1(nlat, nlon, abc, abc(iw1), abc(iw2))

    contains

        subroutine rabcp1(nlat, nlon, a, b, c)
            !
            ! Purpose:
            !
            ! Coefficients a, b, and c for computing pbar(m, n, theta) are
            ! stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
            !
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            real (wp),    intent (in out) :: a(*)
            real (wp),    intent (in out) :: b(*)
            real (wp),    intent (in out) :: c(*)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip) :: m, n, ns, mp1, np1, mp3, mmax
            real (wp)    :: cn, fm, fn
            real (wp)    :: tm, tn, fnmm, fnpm, temp
            !----------------------------------------------------------------------

            mmax = min(nlat, nlon/2+1)

            outer_loop: do mp1=3, mmax

                m = mp1-1
                ns = ((m-2)*(nlat+nlat-m-1))/2+1
                fm = real(m, kind=wp)
                tm = 2.0_wp * fm
                temp = tm*(tm-1.0_wp)
                a(ns) = sqrt((tm+1.0_wp)*(tm-2.0_wp)/temp)
                c(ns) = sqrt(2.0_wp/temp)

                if (m == nlat-1) cycle outer_loop

                ns = ns+1
                temp = tm*(tm+1.0_wp)
                a(ns) = sqrt((tm+3.0_wp)*(tm-2.0_wp)/temp)
                c(ns) = sqrt(6.0_wp/temp)
                mp3 = m+3

                if (mp3 > nlat) cycle outer_loop

                do np1=mp3, nlat
                    n = np1-1
                    ns = ns+1
                    fn = real(n, kind=wp)
                    tn = 2.0_wp * fn
                    cn = (tn+1.0_wp)/(tn-3.0_wp)
                    fnpm = fn+fm
                    fnmm = fn-fm
                    temp = fnpm*(fnpm-1.0_wp)
                    a(ns) = sqrt(cn*(fnpm-3.0_wp)*(fnpm-2.0_wp)/temp)
                    b(ns) = sqrt(cn*fnmm*(fnmm-1.0_wp)/temp)
                    c(ns) = sqrt((fnmm+1.0_wp)*(fnmm+2.0_wp)/temp)
                end do
            end do outer_loop

        end subroutine rabcp1
    end subroutine rabcp


    subroutine sea1(nlat, nlon, imid, z, idz, zin, wzfin, dwork)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        integer (ip), intent (in)     :: imid
        real (wp),    intent (in out) :: z(idz, *)
        integer (ip), intent (in)     :: idz
        real (wp),    intent (in out) :: zin(imid, nlat,3)
        real (wp),    intent (in out) :: wzfin(*)
        real (wp),    intent (in out) :: dwork(*)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: i, m, i3, mn, mp1, np1, mmax
        !----------------------------------------------------------------------

        call zfinit(nlat, nlon, wzfin, dwork)

        mmax = min(nlat, nlon/2+1)

        do mp1=1, mmax
            m = mp1-1

            call zfin(0, nlat, nlon, m, zin, i3, wzfin)

            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                z(mn,1:imid) = zin(:,np1, i3)
            end do
        end do

    end subroutine sea1



    subroutine ses1(nlat, nlon, imid, p, pin, walin, dwork)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        integer (ip), intent (in)     :: imid
        real (wp),    intent (in out) :: p(imid,*)
        real (wp),    intent (in out) :: pin(imid,nlat,3)
        real (wp),    intent (in out) :: walin(*)
        real (wp),    intent (in out) :: dwork(*)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: m, i3, mn, mp1, np1, mmax
        !----------------------------------------------------------------------

        call alinit(nlat, nlon, walin, dwork)

        mmax = min(nlat, nlon/2+1)

        do mp1=1, mmax
            m = mp1-1
            call alin(0, nlat, nlon, m, pin, i3, walin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                p(:, mn) = pin(:, np1, i3)
            end do
        end do

    end subroutine ses1



    subroutine zvinit(nlat, nlon, wzvin, dwork)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        real (wp),    intent (in out) :: wzvin(*)
        real (wp),    intent (in out) :: dwork(nlat+2)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: imid, iw1, iw2
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     the length of wzvin is
        !         2*nlat*imid +3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        !     the length of dwork is nlat+2
        !
        call zvini1(nlat, nlon, imid, wzvin, wzvin(iw1), dwork, dwork(iw2))

    contains

        subroutine zvini1(nlat, nlon, imid, zv, abc, czv, work)
            !
            !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            !     locations where mmax = min(nlat, (nlon+1)/2)
            !     czv and work must each have nlat/2+1  locations
            !
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            integer (ip), intent (in)     :: imid
            real (wp),    intent (in out) :: zv(imid,nlat,2)
            real (wp),    intent (in out) :: abc(*)
            real (wp),    intent (in out) :: czv(nlat/2+1)
            real (wp),    intent (in out) :: work(nlat/2+1)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip)         :: i,m, mdo, mp1, n, np1
            real (wp)            :: dt, th, zvh
            !----------------------------------------------------------------------

            dt = PI/(nlat-1)
            mdo = min(2, nlat, (nlon+1)/2)
            do mp1=1, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dzvk(nlat, m, n, czv, work)
                    do i=1, imid
                        th = real(i-1, kind=wp)*dt
                        call dzvt(nlat, m, n, th, czv, zvh)
                        zv(i, np1, mp1) = zvh
                    end do
                    zv(1, np1, mp1) = 0.5_wp*zv(1, np1, mp1)
                end do
            end do

            call rabcv(nlat, nlon, abc)

        end subroutine zvini1
    end subroutine zvinit



    subroutine zwinit(nlat, nlon, wzwin, dwork)
        !
        ! Purpose:
        !
        ! The length of wzvin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        ! The length of dwork is nlat+2
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        real (wp),    intent (in out) :: wzwin(2*nlat*((nlat+1)/2)+3*((nlat-3)*nlat+2)/2)
        real (wp),    intent (in out) :: dwork(nlat+2)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: imid, iw1, iw2
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2

        call zwini1(nlat, nlon, imid, wzwin, wzwin(iw1), dwork, dwork(iw2))

    contains

        subroutine zwini1(nlat, nlon, imid, zw, abc, czw, work)
            !
            !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            !     locations where mmax = min(nlat, (nlon+1)/2)
            !     czw and work must each have nlat+1 locations
            !
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: nlon
            integer (ip), intent (in)     :: imid
            real (wp),    intent (in out) :: zw(imid,nlat,2)
            real (wp),    intent (in out) :: abc(*)
            real (wp),    intent (in out) :: czw(nlat+1)
            real (wp),    intent (in out) :: work(nlat+1)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip)         :: i, m, mdo, mp1, n, np1
            real (wp)            :: dt, th, zwh
            !----------------------------------------------------------------------

            dt = pi/(nlat-1)
            mdo = min(3, nlat, (nlon+1)/2)

            if (mdo < 2) return

            do mp1=2, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dzwk(nlat, m, n, czw, work)
                    do i=1, imid
                        th = real(i-1, kind=wp)*dt
                        call dzwt(nlat, m, n, th, czw, zwh)
                        zw(i, np1, m) = zwh
                    end do
                    zw(1, np1, m) = 0.5_wp*zw(1, np1, m)
                end do
            end do

            call rabcw(nlat, nlon, abc)

        end subroutine zwini1

    end subroutine zwinit



    subroutine zvin(ityp, nlat, nlon, m, zv, i3, wzvin)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: ityp
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip), intent (in)  :: m
        real (wp),    intent (out) :: zv(*)
        integer (ip), intent (out) :: i3
        real (wp),    intent (in)  :: wzvin(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid
        integer (ip) :: iw1, iw2, iw3, iw4
        integer (ip) :: labc, lim, mmax
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon+1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wzvin is 2*lim+3*labc
        !
        call zvin1(ityp, nlat, m, zv, imid, i3, wzvin, wzvin(iw1), wzvin(iw2), &
            wzvin(iw3), wzvin(iw4))

    contains

        subroutine zvin1(ityp, nlat, m, zv, imid, i3, zvz, zv1, a, b, c)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ityp
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: m
            real (wp),    intent (out)    :: zv(imid, nlat, 3)
            integer (ip), intent (in)     :: imid
            integer (ip), intent (in out) :: i3
            real (wp),    intent (in)     :: zvz(imid, *)
            real (wp),    intent (in)     :: zv1(imid, *)
            real (wp),    intent (in)     :: a(*)
            real (wp),    intent (in)     :: b(*)
            real (wp),    intent (in)     :: c(*)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip)       :: i, ihold
            integer (ip)       :: np1, ns, nstp, nstrt
            integer (ip), save :: i1, i2
            !----------------------------------------------------------------------

            ihold = i1
            i1 = i2
            i2 = i3
            i3 = ihold
            if (m < 1) then
                i1 = 1
                i2 = 2
                i3 = 3
                do np1=1, nlat
                    do i=1, imid
                        zv(i, np1, i3) = zvz(i, np1)
                    end do
                end do
            else if (m == 1) then
                do np1=2, nlat
                    do i=1, imid
                        zv(i, np1, i3) = zv1(i, np1)
                    end do
                end do
            else
                ns = ((m-2)*(nlat+nlat-m-1))/2+1

                if (ityp /= 1) then
                    do i=1, imid
                        zv(i, m+1, i3) = a(ns)*zv(i, m-1, i1)-c(ns)*zv(i, m+1, i1)
                    end do
                end if

                if (m == nlat-1) return

                if (ityp /= 2) then
                    ns = ns+1
                    do i=1, imid
                        zv(i, m+2, i3) = a(ns)*zv(i, m, i1)-c(ns)*zv(i, m+2, i1)
                    end do
                end if

                nstrt = m+3

                if (ityp == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2

                if (ityp == 0) nstp = 1

                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    do i=1, imid
                        zv(i, np1, i3) = a(ns)*zv(i, np1-2, i1)+b(ns)*zv(i, np1-2, i3) &
                            -c(ns)*zv(i, np1, i1)
                    end do
                end do
            end if

        end subroutine zvin1

    end subroutine zvin


    subroutine zwin(ityp, nlat, nlon, m, zw, i3, wzwin)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: ityp
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        integer (ip), intent (in)     :: m
        real (wp),    intent (out)    :: zw(*)
        integer (ip), intent (in out) :: i3
        real (wp),    intent (in)     :: wzwin(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid, iw1, iw2, iw3, iw4
        integer (ip) :: labc, lim, mmax
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon+1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wzwin is 2*lim+3*labc
        !
        call zwin1(ityp, nlat, m, zw, imid, i3, wzwin, wzwin(iw1), wzwin(iw2), &
            wzwin(iw3), wzwin(iw4))

    contains

        subroutine zwin1(ityp, nlat, m, zw, imid, i3, zw1, zw2, a, b, c)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ityp
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: m
            real (wp),    intent (out)    :: zw(imid, nlat, 3)
            integer (ip), intent (in)     :: imid
            integer (ip), intent (in out) :: i3
            real (wp),    intent (in)     :: zw1(imid,*)
            real (wp),    intent (in)     :: zw2(imid,*)
            real (wp),    intent (in)     :: a(*)
            real (wp),    intent (in)     :: b(*)
            real (wp),    intent (in)     :: c(*)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip)       :: i, ihold
            integer (ip)       :: np1, ns, nstp, nstrt
            integer (ip), save :: i1, i2
            !----------------------------------------------------------------------

            ihold = i1
            i1 = i2
            i2 = i3
            i3 = ihold
            if (m < 2) then
                i1 = 1
                i2 = 2
                i3 = 3
                do np1=2, nlat
                    do i=1, imid
                        zw(i, np1, i3) = zw1(i, np1)
                    end do
                end do
            else if (m == 2) then
                do np1=3, nlat
                    do i=1, imid
                        zw(i, np1, i3) = zw2(i, np1)
                    end do
                end do
            else
                ns = ((m-2)*(2*nlat-m-1))/2+1

                if (ityp /= 1) then
                    do i=1, imid
                        zw(i, m+1, i3) = a(ns)*zw(i, m-1, i1)-c(ns)*zw(i, m+1, i1)
                    end do
                end if

                if (m == nlat-1) return

                if (ityp /= 2) then
                    ns = ns+1
                    do i=1, imid
                        zw(i, m+2, i3) = a(ns)*zw(i, m, i1)-c(ns)*zw(i, m+2, i1)
                    end do
                end if

                nstrt = m+3

                if (ityp == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2

                if (ityp == 0) nstp = 1

                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    do i=1, imid
                        zw(i, np1, i3) = a(ns)*zw(i, np1-2, i1)+b(ns)*zw(i, np1-2, i3) &
                            -c(ns)*zw(i, np1, i1)
                    end do
                end do
            end if

        end subroutine zwin1
    end subroutine zwin

    subroutine vbinit(nlat, nlon, wvbin, dwork)
        !
        ! Remark:
        !
        ! The length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        ! The length of dwork is nlat+2
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        real (wp),    intent (out) :: wvbin(2*nlat*((nlat+1)/2)+3*((nlat-3)*nlat+2)/2)
        real (wp),    intent (out) :: dwork(nlat+2)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid, iw1, iw2
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2 + 2

        call vbini1(nlat, nlon, imid, wvbin, wvbin(iw1), dwork, dwork(iw2))

    contains

        subroutine vbini1(nlat, nlon, imid, vb, abc, cvb, work)
            !
            ! Remarks:
            !
            !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            !     locations where mmax = min(nlat, (nlon+1)/2)
            !     cvb and work must each have nlat+1 locations
            !
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)  :: nlat
            integer (ip), intent (in)  :: nlon
            integer (ip), intent (in)  :: imid
            real (wp),    intent (out) :: vb(imid, nlat, 2)
            real (wp),    intent (out) :: abc(*)
            real (wp),    intent (out) :: cvb(nlat+1)
            real (wp),    intent (out) :: work(nlat+1)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip)    :: i, m, mdo, mp1, n, np1
            real (wp)       :: dth, theta, vbh
            !----------------------------------------------------------------------

            dth = pi/(nlat-1)
            mdo = min(2, nlat, (nlon+1)/2)
            do mp1=1, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dvbk(m, n, cvb, work)
                    do  i=1, imid
                        theta = real(i-1, kind=wp)*dth
                        call dvbt(m, n, theta, cvb, vbh)
                        vb(i, np1, mp1) = vbh
                    end do
                end do
            end do

            call rabcv(nlat, nlon, abc)

        end subroutine vbini1

    end subroutine vbinit


    subroutine wbinit (nlat, nlon, wwbin, dwork)
        !
        ! Remark:
        !
        ! The length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        ! The length of dwork is nlat+2
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        real (wp),    intent (out) :: wwbin(2*nlat*((nlat+1)/2)+3*((nlat-3)*nlat+2)/2)
        real (wp),    intent (out) :: dwork(nlat+2)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid, iw1, iw2
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2 + 2

        call wbini1(nlat, nlon, imid, wwbin, wwbin(iw1), dwork, dwork(iw2))

    contains

        subroutine wbini1(nlat, nlon, imid, wb, abc, cwb, work)
            !
            ! Remarks:
            !
            ! abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            ! locations where mmax = min(nlat, (nlon+1)/2)
            ! cwb and work must each have nlat/2+1 locations
            !
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)  :: nlat
            integer (ip), intent (in)  :: nlon
            integer (ip), intent (in)  :: imid
            real (wp),    intent (out) :: wb(imid, nlat, 2)
            real (wp),    intent (out) :: abc(*)
            real (wp),    intent (out) :: cwb(nlat/2+1)
            real (wp),    intent (out) :: work(nlat/2+1)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip)         :: i, m, mdo, mp1, n, np1
            real (wp)            :: dth, wbh, theta
            !----------------------------------------------------------------------

            dth = pi/(nlat-1)
            mdo = min(3, nlat, (nlon+1)/2)

            if (mdo >= 2) then
                do mp1=2, mdo
                    m = mp1-1
                    do np1=mp1, nlat
                        n = np1-1
                        call dwbk(m, n, cwb, work)
                        do i=1, imid
                            theta = real(i-1, kind=wp)*dth
                            call dwbt(m, n, theta, cwb, wbh)
                            wb(i, np1, m) = wbh
                        end do
                    end do
                end do

                call rabcw(nlat, nlon, abc)
            end if

        end subroutine wbini1

    end subroutine wbinit



    subroutine vbin(ityp, nlat, nlon, m, vb, i3, wvbin)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: ityp
        integer (ip), intent (in)     :: nlat
        integer (ip), intent (in)     :: nlon
        integer (ip), intent (in)     :: m
        real (wp),    intent (out)    :: vb(*)
        integer (ip), intent (in out) :: i3
        real (wp),    intent (in)     :: wvbin(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid
        integer (ip) :: iw1, iw2, iw3, iw4
        integer (ip) :: labc, lim, mmax
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon+1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wvbin is 2*lim+3*labc
        !
        call vbin1(ityp, nlat, m, vb, imid, i3, wvbin, wvbin(iw1), wvbin(iw2), &
            wvbin(iw3), wvbin(iw4))

    contains

        subroutine vbin1(ityp, nlat, m, vb, imid, i3, vbz, vb1, a, b, c)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ityp
            integer (ip), intent (in)     :: nlat
            integer (ip), intent (in)     :: m
            real (wp),    intent (out)    :: vb(imid, nlat, 3)
            integer (ip), intent (in)     :: imid
            integer (ip), intent (in out) :: i3
            real (wp),    intent (in)     :: vbz(imid,*)
            real (wp),    intent (in)     :: vb1(imid,*)
            real (wp),    intent (in)     :: a(*)
            real (wp),    intent (in)     :: b(*)
            real (wp),    intent (in)     :: c(*)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip)       :: i, ihold
            integer (ip)       :: np1, ns, nstp, nstrt
            integer (ip), save :: i1, i2
            !----------------------------------------------------------------------

            ihold = i1
            i1 = i2
            i2 = i3
            i3 = ihold
            if (m < 1) then
                i1 = 1
                i2 = 2
                i3 = 3
                do np1=1, nlat
                    do i=1, imid
                        vb(i, np1, i3) = vbz(i, np1)
                    end do
                end do
            else if (m == 1) then
                do np1=2, nlat
                    do i=1, imid
                        vb(i, np1, i3) = vb1(i, np1)
                    end do
                end do
            else
                ns = ((m-2)*(nlat+nlat-m-1))/2+1

                if (ityp /= 1) then
                    do i=1, imid
                        vb(i, m+1, i3) = a(ns)*vb(i, m-1, i1)-c(ns)*vb(i, m+1, i1)
                    end do
                end if

                if (m == nlat-1) return

                if (ityp /= 2) then
                    ns = ns+1
                    do i=1, imid
                        vb(i, m+2, i3) = a(ns)*vb(i, m, i1)-c(ns)*vb(i, m+2, i1)
                    end do
                end if

                nstrt = m+3

                if (ityp == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2

                if (ityp == 0) nstp = 1

                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    do i=1, imid
                        vb(i, np1, i3) = a(ns)*vb(i, np1-2, i1)+b(ns)*vb(i, np1-2, i3) &
                            -c(ns)*vb(i, np1, i1)
                    end do
                end do
            end if

        end subroutine vbin1

    end subroutine vbin


    subroutine wbin(ityp, nlat, nlon, m, wb, i3, wwbin)

        integer (ip) :: i3
        integer (ip) :: imid
        integer (ip), intent (in) :: ityp
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip) :: iw3
        integer (ip) :: iw4
        integer (ip) :: labc
        integer (ip) :: lim
        integer (ip) :: m
        integer (ip) :: mmax
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        real (wp) :: wb(*)
        real (wp) :: wwbin(*)


        imid = (nlat+1)/2
        lim = nlat*imid
        mmax = min(nlat, (nlon+1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = lim+1
        iw2 = iw1+lim
        iw3 = iw2+labc
        iw4 = iw3+labc
        !
        !     the length of wwbin is 2*lim+3*labc
        !
        call wbin1(ityp, nlat, m, wb, imid, i3, wwbin, wwbin(iw1), wwbin(iw2), &
            wwbin(iw3), wwbin(iw4))

    contains

        subroutine wbin1(ityp, nlat, m, wb, imid, i3, wb1, wb2, a, b, c)

            real (wp) :: a(*)
            real (wp) :: b(*)
            real (wp) :: c(*)
            integer (ip) :: i
            integer (ip) :: i3
            integer (ip) :: ihold
            integer (ip) :: imid
            integer (ip), intent (in) :: ityp
            integer (ip) :: m
            integer (ip), intent (in) :: nlat
            integer (ip) :: np1
            integer (ip) :: ns
            integer (ip) :: nstp
            integer (ip) :: nstrt
            real (wp) :: wb(imid, nlat, 3)
            real (wp) :: wb1(imid, *)
            real (wp) ::  wb2(imid, *)
            integer (ip), save :: i1, i2

            ihold = i1
            i1 = i2
            i2 = i3
            i3 = ihold
            if (m < 2) then
                i1 = 1
                i2 = 2
                i3 = 3
                do np1=2, nlat
                    do i=1, imid
                        wb(i, np1, i3) = wb1(i, np1)
                    end do
                end do
            else if (m == 2) then
                do np1=3, nlat
                    do i=1, imid
                        wb(i, np1, i3) = wb2(i, np1)
                    end do
                end do
            else
                ns = ((m-2)*(nlat+nlat-m-1))/2+1

                if (ityp /= 1) then
                    do i=1, imid
                        wb(i, m+1, i3) = a(ns)*wb(i, m-1, i1)-c(ns)*wb(i, m+1, i1)
                    end do
                end if

                if (m == nlat-1) return

                if (ityp /= 2) then
                    ns = ns+1
                    do i=1, imid
                        wb(i, m+2, i3) = a(ns)*wb(i, m, i1)-c(ns)*wb(i, m+2, i1)
                    end do
                end if

                nstrt = m+3

                if (ityp == 1) nstrt = m+4

                if (nstrt > nlat) return

                nstp = 2
                if (ityp == 0) nstp = 1
                do np1=nstrt, nlat, nstp
                    ns = ns+nstp
                    do i=1, imid
                        wb(i, np1, i3) = a(ns)*wb(i, np1-2, i1)+b(ns)*wb(i, np1-2, i3) &
                            -c(ns)*wb(i, np1, i1)
                    end do
                end do
            end if


        end subroutine wbin1

    end subroutine wbin



    subroutine dzvk(nlat, m, n, czv, work)

        integer (ip) :: i
        integer (ip) :: id
        integer (ip) :: k
        integer (ip) :: kdo
        integer (ip) :: lc
        integer (ip) :: m
        integer (ip) :: n
        integer (ip), intent (in) :: nlat
        real (wp) :: work(nlat/2 + 1)
        real (wp) :: czv(*)
        real (wp) :: sc1, summation, t1, t2
        !
        !     subroutine dzvk computes the coefficients in the trigonometric
        !     expansion of the quadrature function zvbar(n, m, theta)
        !
        !     input parameters
        !
        !     nlat      the number of colatitudes including the poles.
        !
        !     n      the degree (subscript) of wbarv(n, m, theta)
        !
        !     m      the order (superscript) of wbarv(n, m, theta)
        !
        !     work   a work array with at least nlat/2+1 locations
        !
        !     output parameter
        !
        !     czv     the fourier coefficients of zvbar(n, m, theta).
        !

        if (n <= 0) return

        lc = (nlat+1)/2
        sc1 = 2.0/(nlat-1)

        call dvbk(m, n, work, czv)

        select case (mod(n, 2))
            case (0)
                select case (mod(m, 2))
                    case (0)
                        !
                        !==> n even, m even
                        !
                        kdo = n/2
                        do id=1, lc
                            i = id+id-2
                            summation = 0.0
                            do k=1, kdo
                                t1 = 1.0-(k+k+i)**2
                                t2 = 1.0-(k+k-i)**2
                                summation = summation+work(k)*(t1-t2)/(t1*t2)
                            end do
                            czv(id) = sc1*summation
                        end do
                    case (1)
                        !
                        !==> n even, m odd
                        !
                        kdo = n/2
                        do id=1, lc
                            i = 2*id-2
                            summation = 0.0
                            do k=1, kdo
                                t1 = 1.0-(k+k+i)**2
                                t2 = 1.0-(k+k-i)**2
                                summation = summation+work(k)*(t1+t2)/(t1*t2)
                            end do
                            czv(id) = sc1*summation
                        end do
                        return
                end select
            case default
                select case (mod(m, 2))
                    case (0)
                        !
                        !==> n odd, m even
                        !
                        kdo = (n+1)/2
                        do id=1, lc
                            i = 2*id-3
                            summation = 0.0
                            do k=1, kdo
                                t1 = 1.0-(k+k-1+i)**2
                                t2 = 1.0-(k+k-1-i)**2
                                summation = summation+work(k)*(t1-t2)/(t1*t2)
                            end do
                            czv(id) = sc1*summation
                        end do
                    case (1)
                        !
                        !==> n odd, m odd
                        !
                        kdo = (n+1)/2
                        do id=1, lc
                            i = 2*id-1
                            summation = 0.0
                            do k=1, kdo
                                t1 = 1.0-(k+k-1+i)**2
                                t2 = 1.0-(k+k-1-i)**2
                                summation = summation+work(k)*(t1+t2)/(t1*t2)
                            end do
                            czv(id) = sc1*summation
                        end do
                end select
        end select

    end subroutine dzvk



    subroutine dzvt(nlat, m, n, th, czv, zvh)

        integer (ip) :: k
        integer (ip) :: lc
        integer (ip) :: lq
        integer (ip) :: ls
        integer (ip) :: m
        integer (ip) :: n
        integer (ip), intent (in) :: nlat
        real (wp) :: czv(*)
        real (wp) :: th, zvh, cost, sint, cdt, sdt, temp
        !
        !     subroutine dzvt tabulates the function zvbar(n, m, theta)
        !     at theta = th in real
        !
        !     input parameters
        !
        !     nlat      the number of colatitudes including the poles.
        !
        !     n      the degree (subscript) of zvbar(n, m, theta)
        !
        !     m      the order (superscript) of zvbar(n, m, theta)
        !
        !     czv     the fourier coefficients of zvbar(n, m, theta)
        !             as computed by subroutine zwk.
        !
        !     output parameter
        !
        !     zvh     zvbar(m, n, theta) evaluated at theta = th
        !

        zvh = 0.0

        if (n <= 0) return

        lc = (nlat+1)/2
        lq = lc-1
        ls = lc-2
        cost = cos(th)
        sint = sin(th)
        cdt = cost**2-sint**2
        sdt = 2.0*sint*cost

        select case (mod(nlat, 2))
            case(0) ! nlat even
                select case (mod(n, 2))
                    case (0) ! n even
                        cost = cdt
                        sint = sdt
                        select case (mod(m, 2))
                            case (0) ! m even
                                !
                                !==> nlat even n even  m even
                                !
                                do k=1, lq
                                    zvh = zvh+czv(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                !
                                !==> nlat even n even m odd
                                !
                                zvh = 0.5*czv(1)
                                do k=2, lc
                                    zvh = zvh+czv(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                        end select
                    case (1) ! n odd
                        select case (mod(m, 2))
                            case (0) ! m even
                                !
                                !==> nlat even n odd  m even
                                !
                                do k=1, lq
                                    zvh = zvh+czv(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                !
                                !==> nlat even  n odd  m odd
                                !
                                zvh = 0.5*czv(lc)*cos(real(nlat-1)*th)

                                do k=1, lq
                                    zvh = zvh+czv(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                        end select
                end select
            case (1) ! nlat odd
                select case (mod(n, 2))
                    case (0) ! n even
                        cost = cdt
                        sint = sdt
                        select case (mod(m, 2))
                            case (0) ! m even
                                !
                                !==> nlat odd  n even  m even
                                !
                                do k=1, ls
                                    zvh = zvh+czv(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                zvh = 0.5*czv(1)
                                do k=2, lq
                                    zvh = zvh+czv(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                                zvh = zvh+0.5*czv(lc)*cos((nlat-1)*th)
                        end select
                    case (1) ! n odd
                        select case (mod(m, 2))
                            case (0) ! m even
                                !
                                !==> nlat odd n odd m even
                                !
                                do k=1, lq
                                    zvh = zvh+czv(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                !
                                !==> nlat odd n odd m odd
                                !
                                do k=1, lq
                                    zvh = zvh+czv(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                        end select
                end select
        end select

    end subroutine dzvt


    subroutine dzwk(nlat, m, n, czw, work)

        integer (ip) :: i
        integer (ip) :: id
        integer (ip) :: k
        integer (ip) :: kdo
        integer (ip) :: kp1
        integer (ip) :: lc
        integer (ip) :: m
        integer (ip) :: n
        integer (ip), intent (in) :: nlat
        real (wp) :: czw(*)
        real (wp) :: work(nlat/2+1)
        real (wp) :: sc1, summation, t1, t2

        !
        !     subroutine dzwk computes the coefficients in the trigonometric
        !     expansion of the quadrature function zwbar(n, m, theta)
        !
        !     input parameters
        !
        !     nlat      the number of colatitudes including the poles.0
        !
        !     n      the degree (subscript) of zwbar(n, m, theta)
        !
        !     m      the order (superscript) of zwbar(n, m, theta)
        !
        !     work   a work array with at least nlat/2+1 locations
        !
        !     output parameter
        !
        !     czw     the fourier coefficients of zwbar(n, m, theta).0
        !

        if (n <= 0) return

        lc = (nlat+1)/2
        sc1 = 2.0/(nlat-1)

        call dwbk(m, n, work, czw)

        select case (mod(n,2))
            case (0) ! n even
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even, m even
                        !
                        kdo = n/2
                        do id=1, lc
                            i = 2*id-3
                            summation = 0.0
                            do k=1, kdo
                                t1 = 1.0-(k+k-1+i)**2
                                t2 = 1.0-(k+k-1-i)**2
                                summation = summation+work(k)*(t1-t2)/(t1*t2)
                            end do
                            czw(id) = sc1*summation
                        end do
                    case (1) ! m odd
                        !
                        !==> n even, m odd
                        !
                        kdo = n/2
                        do id=1, lc
                            i = 2*id-1
                            summation = 0.0
                            do k=1, kdo
                                t1 = 1.0-(k+k-1+i)**2
                                t2 = 1.0-(k+k-1-i)**2
                                summation = summation+work(k)*(t1+t2)/(t1*t2)
                            end do
                            czw(id) = sc1*summation
                        end do
                end select
            case (1) ! n odd
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n odd, m even
                        !
                        kdo = (n-1)/2
                        do id=1, lc
                            i = 2*id-2
                            summation = 0.0
                            do k=1, kdo
                                t1 = 1.0-(k+k+i)**2
                                t2 = 1.0-(k+k-i)**2
                                summation = summation+work(k)*(t1-t2)/(t1*t2)
                            end do
                            czw(id) = sc1*summation
                        end do
                    case (1) ! m odd
                        !
                        !==> n odd, m odd
                        !
                        kdo = (n+1)/2
                        do id=1, lc
                            i = 2*id-2
                            summation = work(1)/(1.0-i**2)

                            if (kdo >= 2) then
                                do kp1=2, kdo
                                    k = kp1-1
                                    t1 = 1.0-(2*k+i)**2
                                    t2 = 1.0-(2*k-i)**2
                                    summation = summation+work(kp1)*(t1+t2)/(t1*t2)
                                end do
                            end if

                            czw(id) = sc1*summation
                        end do
                end select
        end select

    end subroutine dzwk


    subroutine dzwt(nlat, m, n, th, czw, zwh)

        integer (ip) :: k
        integer (ip) :: lc
        integer (ip) :: lq
        integer (ip) :: ls
        integer (ip) :: m
        integer (ip) :: n
        integer (ip), intent (in) :: nlat
        real (wp) :: czw(*)
        real (wp) :: zwh, th, cost, sint, cdt, sdt, temp
        !
        !     subroutine dzwt tabulates the function zwbar(n, m, theta)
        !     at theta = th in real
        !
        !     input parameters
        !
        !     nlat      the number of colatitudes including the poles.
        !            nlat must be an odd integer
        !
        !     n      the degree (subscript) of zwbar(n, m, theta)
        !
        !     m      the order (superscript) of zwbar(n, m, theta)
        !
        !     czw     the fourier coefficients of zwbar(n, m, theta)
        !             as computed by subroutine zwk.
        !
        !     output parameter
        !
        !     zwh     zwbar(m, n, theta) evaluated at theta = th
        !

        zwh = 0.0

        if (n <= 0) return

        lc = (nlat+1)/2
        lq = lc-1
        ls = lc-2
        cost = cos(th)
        sint = sin(th)
        cdt = cost**2-sint**2
        sdt = 2.0*sint*cost

        select case (mod(nlat,2))
            case (0) ! nlat even
                select case (mod(n,2))
                    case (0) ! n even
                        select case (mod(m,2))
                            case (0) ! m even
                                !
                                !==> nlat even  n even  m even
                                !
                                do k=1, lq
                                    zwh = zwh+czw(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                   !
                                   !     nlat even  n even  m odd
                                   !
                                zwh = 0.5*czw(lc)*cos(real(nlat-1)*th)
                                do k=1, lq
                                    zwh = zwh+czw(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                        end select
                    case (1) ! n odd
                        cost = cdt
                        sint = sdt
                        select case (mod(m,2))
                            case (0) ! m even
                                !
                                !==> nlat even  n odd  m even
                                !
                                do k=1, lq
                                    zwh = zwh+czw(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                   !
                                   !==> nlat even  n odd  m odd
                                   !
                                zwh = 0.5*czw(1)
                                do k=2, lc
                                    zwh = zwh+czw(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                        end select
                end select
            case (1) ! nlat odd
                select case (mod(n,2))
                    case (0) ! n even
                        select case (mod(m,2))
                            case (0) ! m even
                                !
                                !==> nlat odd  n even  m even
                                !
                                do k=1, lq
                                    zwh = zwh+czw(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                !
                                !==> nlat odd  n even  m odd
                                !
                                do k=1, lq
                                    zwh = zwh+czw(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                        end select
                    case (1) ! n odd
                        cost = cdt
                        sint = sdt
                        select case (mod(m,2))
                            case (0) ! m even
                                !
                                !==> nlat odd  n odd  m even
                                !
                                do k=1, ls
                                    zwh = zwh+czw(k+1)*sint
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                            case (1) ! m odd
                                 !
                                 !==> nlat odd  n odd  m odd
                                 !
                                zwh = 0.5*czw(1)
                                do k=2, lq
                                    zwh = zwh+czw(k)*cost
                                    temp = cdt*cost-sdt*sint
                                    sint = sdt*cost+cdt*sint
                                    cost = temp
                                end do
                                zwh = zwh+0.5*czw(lc)*cos(real(nlat-1)*th)
                        end select
                end select
        end select

    end subroutine dzwt



    subroutine dvbk(m, n, cv, work)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: cv(*)
        real (wp),    intent (out) :: work(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: l, ncv
        real (wp)    :: srnp1, fn, fk, cf
        !----------------------------------------------------------------------

        cv(1) = 0.0_wp

        if (n <= 0) return

        fn = n
        srnp1 = sqrt(fn * (fn + 1.0_wp))
        cf = 2.0_wp*real(m, kind=wp)/srnp1

        call dnlfk(m, n, work)

        select case (mod(n,2))
            case (0) ! n even
                ncv = n/2
                if (ncv == 0) return
                fk = 0.0
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even m even
                        !
                        do l=1, ncv
                            fk = fk+2.0
                            cv(l) = -fk*work(l+1)/srnp1
                        end do
                    case (1) ! m odd
                        !
                        !==> n even m odd
                        !
                        do l=1, ncv
                            fk = fk+2.0_wp
                            cv(l) = fk*work(l)/srnp1
                        end do
                end select
            case (1) ! n odd
                ncv = (n+1)/2
                fk = -1.0_wp
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !     n odd m even
                        !
                        do l=1, ncv
                            fk = fk+2.0_wp
                            cv(l) = -fk*work(l)/srnp1
                        end do
                    case (1) ! m odd
                        !
                        !==>     n odd m odd
                        !
                        do l=1, ncv
                            fk = fk+2.0_wp
                            cv(l) = fk*work(l)/srnp1
                        end do
                end select
        end select

    end subroutine dvbk


    subroutine dwbk(m, n, cw, work)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: cw(*)
        real (wp),    intent (out) :: work(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: l
        real (wp)    :: fn, cf, srnp1
        !----------------------------------------------------------------------

        cw(1) = 0.0_wp

        if (n <= 0 .or. m <= 0) return

        fn = n
        srnp1 = sqrt(fn * (fn + 1.0_wp))
        cf = 2.0_wp * real(m, kind=wp)/srnp1

        call dnlfk(m, n, work)

        if (m == 0) return

        select case (mod(n,2))
            case (0) ! n even
                l = n/2
                if (l == 0) return
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even m even
                        !
                        cw(l) = -cf*work(l+1)
                        do
                            l = l-1
                            if (l <= 0) exit
                            cw(l) = cw(l+1)-cf*work(l+1)
                        end do
                    case (1) ! m odd
                             !
                             !     n even m odd
                             !
                        cw(l) = cf*work(l)
                        do
                            l = l-1
                            if (l <= 0) exit
                            cw(l) = cw(l+1)+cf*work(l)
                        end do
                end select
            case (1) ! n odd
                select case (mod(m,2))
                    case (0) ! m even
                        l = (n-1)/2
                        if (l == 0) return
                        !
                        !==> n odd m even
                        !
                        cw(l) = -cf*work(l+1)
                        do
                            l = l-1
                            if (l <= 0) exit
                            cw(l) = cw(l+1)-cf*work(l+1)
                        end do
                    case (1) ! m odd
                        !
                        !==> n odd m odd
                        !
                        l = (n+1)/2
                        cw(l) = cf*work(l)
                        do
                            l = l-1
                            if (l <= 0) exit
                            cw(l) = cw(l+1)+cf*work(l)
                        end do
                end select
        end select

    end subroutine dwbk


    subroutine dvbt(m, n, theta, cv, vh)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (in)  :: theta
        real (wp),    intent (out) :: cv(*)
        real (wp),    intent (out) :: vh
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: k, ncv
        real (wp)    :: cost, sint, cdt, sdt, temp
        !----------------------------------------------------------------------

        vh = 0.0_wp

        if (n == 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = cost**2-sint**2
        sdt = 2.0*sint*cost

        select case (mod(n,2))
            case (0) ! n even
                cost = cdt
                sint = sdt
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even m even
                        !
                        ncv = n/2
                        do k=1, ncv
                            vh = vh+cv(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                        !
                        !==> n even  m odd
                        !
                        ncv = n/2
                        do k=1, ncv
                            vh = vh+cv(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                end select
            case (1) ! n odd
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n odd m even
                        !
                        ncv = (n+1)/2
                        do k=1, ncv
                            vh = vh+cv(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                        !
                        !==> n odd m odd
                        !
                        ncv = (n+1)/2
                        do k=1, ncv
                            vh = vh+cv(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                end select
        end select

    end subroutine dvbt



    subroutine dwbt(m, n, theta, cw, wh)

        integer (ip) :: k
        integer (ip) :: m
        integer (ip) :: n
        integer (ip) :: ncw
        real (wp) :: cw(*)
        real (wp) :: theta, wh, cost, sint, cdt, sdt, temp

        wh = 0.0

        if (n <= 0 .or. m <= 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = cost*cost-sint*sint
        sdt = 2.0*sint*cost

        select case (mod(n,2))
            case (0) ! n even
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even  m even
                        !
                        ncw = n/2
                        do k=1, ncw
                            wh = wh+cw(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                         !
                         !==> n even m odd
                         !
                        ncw = n/2
                        do k=1, ncw
                            wh = wh+cw(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                end select
            case (1) ! n odd
                cost = cdt
                sint = sdt
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n odd m even
                        !
                        ncw = (n-1)/2
                        do k=1, ncw
                            wh = wh+cw(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                        !
                        !==> n odd m odd
                        !
                        ncw = (n+1)/2
                        wh = 0.5*cw(1)

                        if (ncw < 2) return

                        do k=2, ncw
                            wh = wh+cw(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do

                end select
        end select

    end subroutine dwbt



    subroutine rabcv(nlat, nlon, abc)

        real (wp) :: abc(*)
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip) :: labc
        integer (ip) :: mmax
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        !
        !     subroutine rabcp computes the coefficients in the recurrence
        !     relation for the functions vbar(m, n, theta). array abc
        !     must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2 locations.
        !
        mmax = min(nlat, (nlon+1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = labc+1
        iw2 = iw1+labc

        call rabcv1(nlat, nlon, abc, abc(iw1), abc(iw2))

    contains

        subroutine rabcv1(nlat, nlon, a, b, c)

            real (wp) :: a(*)
            real (wp) :: b(*)
            real (wp) :: c(*)
            real (wp) :: cn
            real (wp) :: fm
            real (wp) :: fn
            real (wp) :: fnmm
            real (wp) :: fnpm
            integer (ip) :: m
            integer (ip) :: mmax
            integer (ip) :: mp1
            integer (ip) :: mp3
            integer (ip) :: n
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip) :: np1
            integer (ip) :: ns
            real (wp) :: temp
            real (wp) :: tm
            real (wp) :: tn
            real (wp) :: tpn
            !
            !     coefficients a, b, and c for computing vbar(m, n, theta) are
            !     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
            !

            mmax = min(nlat, (nlon+1)/2)

            if (mmax < 3) return

            outer_loop: do mp1=3, mmax
                m = mp1-1
                ns = ((m-2)*(2*nlat-m-1))/2+1
                fm = real(m, kind=wp)
                tm = fm+fm
                temp = tm*(tm-1.0)
                tpn = (fm-2.0)*(fm-1.0)/(fm*(fm+1.0))
                a(ns) = sqrt(tpn*(tm+1.0)*(tm-2.0)/temp)
                c(ns) = sqrt(2.0/temp)
                if (m == nlat-1) cycle outer_loop
                ns = ns+1
                temp = tm*(tm+1.0)
                tpn = (fm-1.0)*fm/((fm+1.0)*(fm+2.0))
                a(ns) = sqrt(tpn*(tm+3.0)*(tm-2.0)/temp)
                c(ns) = sqrt(6.0/temp)
                mp3 = m+3
                if (mp3 > nlat) cycle outer_loop
                do np1=mp3, nlat
                    n = np1-1
                    ns = ns+1
                    fn = real(n)
                    tn = 2.0*fn
                    cn = (tn+1.0)/(tn-3.0)
                    tpn = (fn-2.0)*(fn-1.0)/(fn*(fn + 1.0))
                    fnpm = fn+fm
                    fnmm = fn-fm
                    temp = fnpm*(fnpm-1.0)
                    a(ns) = sqrt(tpn*cn*(fnpm-3.0)*(fnpm-2.0)/temp)
                    b(ns) = sqrt(tpn*cn*fnmm*(fnmm-1.0)/temp)
                    c(ns) = sqrt((fnmm+1.0)*(fnmm+2.0)/temp)
                end do
            end do outer_loop

        end subroutine rabcv1

    end subroutine rabcv


    subroutine rabcw(nlat, nlon, abc)

        real (wp) :: abc(*)
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip) :: labc
        integer (ip) :: mmax
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        !
        !     subroutine rabcw computes the coefficients in the recurrence
        !     relation for the functions wbar(m, n, theta). array abc
        !     must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2 locations.
        !
        mmax = min(nlat, (nlon+1)/2)
        labc = (max(mmax-2, 0)*(2*nlat-mmax-1))/2
        iw1 = labc+1
        iw2 = iw1+labc
        call rabcw1(nlat, nlon, abc, abc(iw1), abc(iw2))

    contains

        subroutine rabcw1(nlat, nlon, a, b, c)

            real (wp) :: a(*)
            real (wp) :: b(*)
            real (wp) :: c(*)
            real (wp) :: cn
            real (wp) :: fm
            real (wp) :: fn
            real (wp) :: fnmm
            real (wp) :: fnpm
            integer (ip) :: m
            integer (ip) :: mmax
            integer (ip) :: mp1
            integer (ip) :: mp3
            integer (ip) :: n
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip) :: np1
            integer (ip) :: ns
            real (wp) :: temp
            real (wp) :: tm
            real (wp) :: tn
            real (wp) :: tph
            real (wp) :: tpn
            !
            !     coefficients a, b, and c for computing wbar(m, n, theta) are
            !     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
            !

            mmax = min(nlat, (nlon+1)/2)

            if (mmax < 4) return

            outer_loop: do mp1=4, mmax
                m = mp1-1
                ns = ((m-2)*(nlat+nlat-m-1))/2+1
                fm = real(m, kind=wp)
                tm = 2.0*fm
                temp = tm*(tm-1.0)
                tpn = (fm-2.0)*(fm-1.0)/(fm*(fm+1.0))
                tph = fm/(fm-2.0)
                a(ns) = tph*sqrt(tpn*(tm+1.0)*(tm-2.0)/temp)
                c(ns) = tph*sqrt(2.0/temp)
                if (m == nlat-1) cycle outer_loop
                ns = ns+1
                temp = tm*(tm+1.0)
                tpn = (fm-1.0)*fm/((fm+1.0)*(fm+2.0))
                tph = fm/(fm-2.0)
                a(ns) = tph*sqrt(tpn*(tm+3.0)*(tm-2.0)/temp)
                c(ns) = tph*sqrt(6.0/temp)
                mp3 = m+3
                if (mp3 > nlat) cycle outer_loop
                do np1=mp3, nlat
                    n = np1-1
                    ns = ns+1
                    fn = real(n)
                    tn = 2.0*fn
                    cn = (tn+1.0)/(tn-3.0)
                    fnpm = fn+fm
                    fnmm = fn-fm
                    temp = fnpm*(fnpm-1.0)
                    tpn = (fn-2.0)*(fn-1.0)/(fn*(fn + 1.0))
                    tph = fm/(fm-2.0)
                    a(ns) = tph*sqrt(tpn*cn*(fnpm-3.0)*(fnpm-2.0)/temp)
                    b(ns) = sqrt(tpn*cn*fnmm*(fnmm-1.0)/temp)
                    c(ns) = tph*sqrt((fnmm+1.0)*(fnmm+2.0)/temp)
                end do
            end do outer_loop

        end subroutine rabcw1

    end subroutine rabcw

    subroutine vtinit (nlat, nlon, wvbin, dwork)

        integer (ip) :: imid
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        real (wp) :: wvbin(*)
        real (wp) :: dwork(nlat+2)


        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of dwork is nlat+2
        !
        call vtini1(nlat, nlon, imid, wvbin, wvbin(iw1), dwork, dwork(iw2))

    contains

        subroutine vtini1(nlat, nlon, imid, vb, abc, cvb, work)
            !
            !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            !     locations where mmax = min(nlat, (nlon+1)/2)
            !     cvb and work must each have nlat/2+1 locations
            !
            real (wp) :: abc(*)
            integer (ip) :: i
            integer (ip) :: imid
            integer (ip) :: m
            integer (ip) :: mdo
            integer (ip) :: mp1
            integer (ip) :: n
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip) :: np1
            real (wp) :: vb(imid, nlat, 2)
            real (wp) :: dt
            real (wp) :: cvb(nlat/2+1)
            real (wp) :: th, vbh
            real (wp) :: work(nlat/2+1)


            dt = PI/(nlat-1)
            mdo = min(2, nlat, (nlon+1)/2)
            do mp1=1, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dvtk(m, n, cvb, work)
                    do i=1, imid
                        th = real(i-1)*dt
                        call dvtt(m, n, th, cvb, vbh)
                        vb(i, np1, mp1) = vbh
                    end do
                end do
            end do
            call rabcv(nlat, nlon, abc)

        end subroutine vtini1

    end subroutine vtinit


    subroutine wtinit(nlat, nlon, wwbin, dwork)

        integer (ip) :: imid
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        real (wp) :: wwbin(*)
        real (wp) :: dwork(nlat+2)

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of dwork is nlat+2
        !
        call wtini1(nlat, nlon, imid, wwbin, wwbin(iw1), dwork, dwork(iw2))

    contains

        subroutine wtini1(nlat, nlon, imid, wb, abc, cwb, work)

            real (wp) :: abc(*)
            integer (ip) :: i
            integer (ip) :: imid
            integer (ip) :: m
            integer (ip) :: mdo
            integer (ip) :: mp1
            integer (ip) :: n
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip) :: np1
            real (wp) :: wb(imid, nlat, 2)
            real (wp) :: dt
            real (wp) :: cwb(nlat/2+1), wbh, th
            real (wp) :: work(nlat/2+1)
            !
            !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            !     locations where mmax = min(nlat, (nlon+1)/2)
            !     cwb and work must each have nlat/2+1 locations
            !

            dt = PI/(nlat-1)
            mdo = min(3, nlat, (nlon+1)/2)
            if (mdo < 2) return
            do mp1=2, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dwtk(m, n, cwb, work)
                    do i=1, imid
                        th = real(i-1)*dt
                        call dwtt(m, n, th, cwb, wbh)
                        wb(i, np1, m) = wbh
                    end do
                end do
            end do

            call rabcw(nlat, nlon, abc)

        end subroutine wtini1

    end subroutine wtinit

    subroutine vtgint (nlat, nlon, theta, wvbin, work)

        integer (ip) :: imid
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        real (wp) :: wvbin(*)
        real (wp) :: theta((nlat+1)/2)
        real (wp) :: work(nlat+2)


        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat+1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call vtgit1(nlat, nlon, imid, theta, wvbin, wvbin(iw1), work, work(iw2))

    contains

        subroutine vtgit1(nlat, nlon, imid, theta, vb, abc, cvb, work)

            real (wp) :: abc(*)
            integer (ip) :: i
            integer (ip) :: imid
            integer (ip) :: m
            integer (ip) :: mdo
            integer (ip) :: mp1
            integer (ip) :: n
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip) :: np1
            real (wp) :: vb(imid, nlat, 2)
            real (wp) :: theta(*), cvb(nlat/2+1), work(nlat/2+1), vbh
            !
            !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            !     locations where mmax = min(nlat, (nlon+1)/2)
            !     cvb and work must each have nlat/2+1   locations
            !

            mdo = min(2, nlat, (nlon+1)/2)
            do mp1=1, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dvtk(m, n, cvb, work)
                    do i=1, imid
                        call dvtt(m, n, theta(i), cvb, vbh)
                        vb(i, np1, mp1) = vbh
                    end do
                end do
            end do

            call rabcv(nlat, nlon, abc)

        end subroutine vtgit1

    end subroutine vtgint

    subroutine wtgint(nlat, nlon, theta, wwbin, work)

        integer (ip) :: imid
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        real (wp) :: wwbin(*)
        real (wp) :: theta((nlat+1)/2)
        real (wp) :: work(nlat+2)


        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat+1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call wtgit1(nlat, nlon, imid, theta, wwbin, wwbin(iw1), work, work(iw2))

    contains

        subroutine wtgit1(nlat, nlon, imid, theta, wb, abc, cwb, work)

            real (wp) :: abc(*)
            integer (ip) :: i
            integer (ip) :: imid
            integer (ip) :: m
            integer (ip) :: mdo
            integer (ip) :: mp1
            integer (ip) :: n
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip) :: np1
            real (wp) :: wb(imid, nlat, 2)
            real (wp) :: theta(*)
            real (wp) :: cwb(nlat/2+1)
            real (wp) :: work(nlat/2+1), wbh
            !
            !     abc must have 3*((nlat-3)*nlat+2)/2 locations
            !     cwb and work must each have nlat/2+1 locations
            !

            mdo = min(3, nlat, (nlon+1)/2)
            if (mdo < 2) return
            do mp1=2, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dwtk(m, n, cwb, work)
                    do i=1, imid
                        call dwtt(m, n, theta(i), cwb, wbh)
                        wb(i, np1, m) = wbh
                    end do
                end do
            end do

            call rabcw(nlat, nlon, abc)

        end subroutine wtgit1

    end subroutine wtgint



    subroutine dvtk(m, n, cv, work)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: cv(*)
        real (wp),    intent (out) :: work(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: l, ncv
        real (wp)    :: fn, fk, cf, srnp1
        !----------------------------------------------------------------------

        cv(1) = 0.0_wp

        if (n <= 0) return

        fn = n
        srnp1 = sqrt(fn * (fn + 1.0_wp))
        cf = 2.0_wp * real(m, kind=wp)/srnp1

        call dnlfk(m, n, work)

        select case (mod(n,2))
            case (0) ! n even
                ncv = n/2
                if (ncv == 0) return
                fk = 0.0
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even m even
                        !
                        do l=1, ncv
                            fk = fk+2.0
                            cv(l) = -(fk**2)*work(l+1)/srnp1
                        end do
                    case (1) ! m odd
                        !
                        !==> n even m odd
                        !
                        do l=1, ncv
                            fk = fk+2.0
                            cv(l) = -(fk**2)*work(l)/srnp1
                        end do
                end select
            case (1) ! n odd
                ncv = (n+1)/2
                fk = -1.0
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n odd m even
                        !
                        do l=1, ncv
                            fk = fk+2.0
                            cv(l) = -(fk**2)*work(l)/srnp1
                        end do
                    case (1) ! m odd
                        !
                        !==> n odd m odd
                        !
                        do l=1, ncv
                            fk = fk+2.0
                            cv(l) = -(fk**2)*work(l)/srnp1
                        end do
                end select
        end select

    end subroutine dvtk



    subroutine dwtk(m, n, cw, work)

        integer (ip) :: l
        integer (ip) :: m
        integer (ip) :: n
        real cw(*), work(*), fn, cf, srnp1

        cw(1) = 0.0

        if (n <= 0 .or. m <= 0) return

        fn = n
        srnp1 = sqrt(fn * (fn + 1.0_wp))
        cf = 2.0_wp * real(m, kind=wp)/srnp1

        call dnlfk(m, n, work)

        if (m == 0) return

        select case (mod(n,2))
            case (0) ! n even
                l = n/2
                if (l == 0) return
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !     n even m even
                        !
                        cw(l) = -cf*work(l+1)
                        do
                            l = l-1
                            if (l <= 0) exit
                            cw(l) = cw(l+1)-cf*work(l+1)
                            cw(l+1) = (2*l+1)*cw(l+1)
                        end do
                    case (1) ! m odd
                        !
                        !==> n even m odd
                        !
                        cw(l) = cf*work(l)
                        do
                            l = l-1
                            if (l < 0) then
                                exit
                            else if (l == 0) then
                                cw(l+1) = -(2*l+1)*cw(l+1)
                            else
                                cw(l) = cw(l+1)+cf*work(l)
                                cw(l+1) = -(2*l+1)*cw(l+1)
                            end if
                        end do
                end select
            case (1) ! n odd
                select case (mod(m,2))
                    case (0) ! m even
                        l = (n-1)/2
                        if (l == 0) return
                        !
                        !==> n odd m even
                        !
                        cw(l) = -cf*work(l+1)
                        do
                            l = l-1
                            if (l < 0) then
                                exit
                            else if (l == 0) then
                                !cw(l) = cw(l+1)-cf*work(l+1)
                                cw(l+1) = (2*l+2)*cw(l+1)
                            else
                                cw(l) = cw(l+1)-cf*work(l+1)
                                cw(l+1) = (2*l+2)*cw(l+1)
                            end if
                        end do
                    case (1) ! m odd
                        !
                        !==> n odd m odd
                        !
                        l = (n+1)/2
                        cw(l) = cf*work(l)
                        do
                            l = l-1
                            if (l < 0) then
                                exit
                            else if (l == 0) then
                                cw(l+1) = -(2*l)*cw(l+1)
                            else
                                cw(l) = cw(l+1)+cf*work(l)
                                cw(l+1) = -(2*l)*cw(l+1)
                            end if
                        end do
                end select
        end select

    end subroutine dwtk



    subroutine dvtt(m, n, theta, cv, vh)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: m
        integer (ip), intent (in)  :: n
        real (wp),    intent (in)  :: theta
        real (wp),    intent (out) :: cv(*)
        real (wp),    intent (out) :: vh
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: k, ncv
        real (wp)    :: cost, sint, cdt, sdt, temp
        !----------------------------------------------------------------------

        vh = 0.0_wp

        if (n == 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = cost**2-sint**2
        sdt = 2.0_wp*sint*cost

        select case (mod(n,2))
            case (0) ! n even
                cost = cdt
                sint = sdt
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even  m even
                        !
                        ncv = n/2
                        do k=1, ncv
                            vh = vh+cv(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                         !
                         !==> n even  m odd
                         !
                        ncv = n/2
                        do k=1, ncv
                            vh = vh+cv(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                end select
            case (1) ! n odd
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !     n odd m even
                        !
                        ncv = (n+1)/2
                        do k=1, ncv
                            vh = vh+cv(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                         !
                         !==> n odd m odd
                         !
                        ncv = (n+1)/2
                        do k=1, ncv
                            vh = vh+cv(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                end select
        end select

    end subroutine dvtt


    subroutine dwtt(m, n, theta, cw, wh)

        integer (ip) :: k
        integer (ip) :: m
        integer (ip) :: n
        integer (ip) :: ncw
        real (wp) :: cw(*)
        real (wp) :: theta, wh, cost, sint, cdt, sdt, temp

        wh = 0.0

        if (n <= 0 .or. m <= 0) return

        cost = cos(theta)
        sint = sin(theta)
        cdt = cost**2-sint**2
        sdt = 2.0*sint*cost

        select case (mod(n,2))
            case (0) ! n even
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n even m even
                        !
                        ncw = n/2
                        do k=1, ncw
                            wh = wh+cw(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                          !
                          !==> n even m odd
                          !
                        ncw = n/2
                        do k=1, ncw
                            wh = wh+cw(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                end select
            case (1) ! n odd
                cost = cdt
                sint = sdt
                select case (mod(m,2))
                    case (0) ! m even
                        !
                        !==> n odd m even
                        !
                        ncw = (n-1)/2
                        do k=1, ncw
                            wh = wh+cw(k)*cost
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                    case (1) ! m odd
                          !
                          !==> n odd m odd
                          !
                        ncw = (n+1)/2
                        wh = 0.0

                        if (ncw < 2) return

                        do k=2, ncw
                            wh = wh+cw(k)*sint
                            temp = cdt*cost-sdt*sint
                            sint = sdt*cost+cdt*sint
                            cost = temp
                        end do
                end select
        end select

    end subroutine dwtt




    subroutine vbgint (nlat, nlon, theta, wvbin, work)

        integer (ip) :: imid
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        real (wp) :: wvbin(*)
        real (wp) :: theta((nlat+1)/2)
        real (wp) :: work(nlat+2)

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat+1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call vbgit1(nlat, nlon, imid, theta, wvbin, wvbin(iw1), work, work(iw2))


    contains

        subroutine vbgit1(nlat, nlon, imid, theta, vb, abc, cvb, work)

            real (wp) :: abc(*)
            integer (ip) :: i
            integer (ip) :: imid
            integer (ip) :: m
            integer (ip) :: mdo
            integer (ip) :: mp1
            integer (ip) :: n
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip) :: np1
            real (wp) :: vb(imid, nlat, 2)
            real (wp) :: cvb(nlat/2+1)
            real (wp) :: theta(*), vbh
            real (wp) :: work(nlat/2+1)
            !
            !     abc must have 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
            !     locations where mmax = min(nlat, (nlon+1)/2)
            !     cvb and work must each have nlat/2+1 locations
            !

            mdo = min(2, nlat, (nlon+1)/2)
            do mp1=1, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dvbk(m, n, cvb, work)
                    do i=1, imid
                        call dvbt(m, n, theta(i), cvb, vbh)
                        vb(i, np1, mp1) = vbh
                    end do
                end do
            end do

            call rabcv(nlat, nlon, abc)

        end subroutine vbgit1

    end subroutine vbgint


    subroutine wbgint(nlat, nlon, theta, wwbin, work)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: nlat
        integer (ip), intent (in) :: nlon
        real (wp),    intent (in) :: theta((nlat+1)/2)
        real (wp)                 :: wwbin(2*nlat*((nlat+1)/2)+3*((nlat-3)*nlat+2)/2)
        real (wp)                 :: work(nlat+2)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip) :: imid, iw1, iw2
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        iw1 = 2*nlat*imid+1
        iw2 = nlat/2+2
        !
        !     theta is a real array with (nlat+1)/2 locations
        !     nlat is the maximum value of n+1
        !     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
        !     the length of work is nlat+2
        !
        call wbgit1(nlat, nlon, imid, theta, wwbin, wwbin(iw1), work, work(iw2))

    contains

        subroutine wbgit1(nlat, nlon, imid, theta, wb, abc, cwb, work)
            !----------------------------------------------------------------------
            ! Dummy arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in) :: nlat
            integer (ip), intent (in) :: nlon
            integer (ip), intent (in) :: imid
            real (wp),    intent (in) :: theta(*)
            real (wp)                 :: wb(imid, nlat, 2)
            real (wp)                 :: abc(3*((nlat-3)*nlat+2)/2)
            real (wp)                 :: cwb(nlat/2+1)
            real (wp)                 :: work(nlat/2+1)
            !----------------------------------------------------------------------
            ! Local variables
            !----------------------------------------------------------------------
            integer (ip) :: i, m, mdo, mp1, n, np1
            real (wp)    :: wbh
            !----------------------------------------------------------------------

            !
            !     abc must have 3*((nlat-3)*nlat+2)/2 locations
            !     cwb and work must each have nlat/2+1 locations
            !
            mdo = min(3, nlat, (nlon+1)/2)
            if (mdo < 2) return
            do mp1=2, mdo
                m = mp1-1
                do np1=mp1, nlat
                    n = np1-1
                    call dwbk(m, n, cwb, work)
                    do i=1, imid
                        call dwbt(m, n, theta(i), cwb, wbh)
                        wb(i, np1, m) = wbh
                    end do
                end do
            end do

            call rabcw(nlat, nlon, abc)

        end subroutine wbgit1

    end subroutine wbgint

end module type_SpherepackAux
