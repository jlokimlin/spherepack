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
!     *       A Package of Fortran77 Subroutines and Programs         *
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
!
!     This version of gaqd implements the method presented in:
!     P. N. swarztrauber, Computing the points and weights for
!     Gauss-Legendre quadrature, SIAM J. Sci. Comput.,
!     24(2002) pp. 945-954.
!
!     The w and lwork arrays are dummy and included only to
!     permit a simple pluggable exchange with the
!     old gaqd in previous versions of SPHEREPACK
!
!
subroutine gaqd(nlat, theta, wts, w, lwork, ierror)
    !
    !     gauss points and weights are computed using the fourier-newton
    !     described in "on computing the points and weights for 
    !     gauss-legendre quadrature", paul n. swarztrauber, siam journal 
    !     on scientific computing that has been accepted for publication.
    !     This routine is faster and more accurate than older program
    !     with the same name.
    !
    !     subroutine gaqd computes the nlat gaussian colatitudes and weights
    !     in real. the colatitudes are in radians and lie in the
    !     in the interval (0, pi).
    !
    !     input parameters
    !
    !     nlat    the number of gaussian colatitudes in the interval (0, pi)
    !             (between the two poles).  nlat must be greater than zero.
    !
    !     w       unused real variable that permits a simple
    !             exchange with the old routine with the same name
    !             in spherepack.
    !
    !     lwork   unused variable that permits a simple exchange with the
    !             old routine with the same name in spherepack.
    !
    !     output parameters
    !
    !     theta   a real array with length nlat
    !             containing the gaussian colatitudes in
    !             increasing radians on the interval (0, pi).
    !
    !     wts     a real array with lenght nlat
    !             containing the gaussian weights.
    !
    !     ierror = 0 no errors
    !            = 1 if nlat <= 0
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)  :: nlat
    real (wp),    intent (out) :: theta(nlat)
    real (wp),    intent (out) :: wts(nlat)
    real (wp),    intent (in)  :: w
    integer (ip), intent (in)  :: lwork
    integer (ip), intent (out) :: ierror
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip)         :: i, idx, it, mnlat, nix
    integer (ip)         :: nhalf, ns2
    real (wp), parameter :: PI = acos(-1.0_wp)
    real (wp), parameter :: HALF_PI = PI/2
    real (wp), parameter :: WP_ZERO = nearest(1.0_wp,1.0_wp)-nearest(1.0_wp,-1.0_wp)
    real (wp)            :: dtheta, dthalf, cmax
    real (wp)            :: zprev, zlast, zero
    real (wp)            :: zhold, pb, dpb, dcor, cz
    real (wp)            :: eps, sgnd
    !----------------------------------------------------------------------

    !
    ! Initialize error flag
    ierror = 0

    ! Check case 1: work space length
    if (nlat <= 0) then
        ierror = 1
        return
    end if

    !
    !==> compute weights and points analytically when nlat=1, 2
    !
    select case (nlat)
        case(1)
            theta(1) = acos(0.0_wp)
            wts(1) = 2.0_wp
            return
        case(2)
            associate( x => sqrt(1.0_wp/3.0_wp) )
                theta(1) = acos(x)
                theta(2) = acos(-x)
            end associate
            wts(1) = 1.0_wp
            wts(2) = 1.0_wp
            return
        case default
            eps = sqrt(epsilon(1.0_wp))
            eps = eps * sqrt(eps)
            mnlat = mod(nlat, 2)
            ns2 = nlat/2
            nhalf = (nlat+1)/2
            idx = ns2+2

            call cpdp(nlat, cz, theta(ns2+1), wts(ns2+1))

            dtheta = HALF_PI/nhalf
            dthalf = dtheta/2
            cmax = 0.2_wp * dtheta
            !
            !==> Estimate first point next to theta = pi/2
            !
            if (mnlat /= 0) then
                zero = HALF_PI-dtheta
                zprev = HALF_PI
                nix = nhalf-1
            else
                zero = HALF_PI-dthalf
                nix = nhalf
            end if

            do while (nix /= 0)
                it = 0
                do while (abs(zero-zlast) > eps*abs(zero) .or. it == 0 )
                    !
                    !==> Increment iterator
                    !
                    it = it+1
                    zlast = zero
                    !
                    !==> newton iterations
                    !
                    call tpdp(nlat, zero, cz, theta(ns2+1), wts(ns2+1), pb, dpb)
                    dcor = pb/dpb
                    sgnd = 1.0_wp

                    if ( dcor /= WP_ZERO ) then
                        sgnd = dcor/abs(dcor)
                    end if

                    dcor = sgnd * min(abs(dcor), cmax)
                    zero = zero-dcor

                end do

                theta(nix) = zero
                zhold = zero
                !    wts(nix) = (nlat+nlat+1)/(dpb*dpb)
                !
                ! ==> yakimiw's formula permits using old pb and dpb
                !
                wts(nix) = (2*nlat+1)/(dpb+pb*cos(zlast)/sin(zlast))**2
                nix = nix-1

                if (nix == nhalf-1)  then
                    zero = 3.0_wp * zero - PI
                end if

                if (nix < nhalf-1)  then
                    zero = 2.0_wp * zero-zprev
                end if

                zprev = zhold

            end do
            !
            !==> Extend points and weights via symmetries
            !
            if (mnlat /= 0) then
                theta(nhalf) = HALF_PI
                call tpdp(nlat, HALF_PI, cz, theta(ns2+1), wts(ns2+1), pb, dpb)
                wts(nhalf) = (2*nlat+1)/(dpb*dpb)
            end if

            do i=1, ns2
                wts(nlat-i+1) = wts(i)
                theta(nlat-i+1) = PI-theta(i)
            end do

            !
            !==> Set weights
            !
            wts = 2.0_wp * wts/sum(wts)
    end select

end subroutine gaqd


pure subroutine cpdp(n, cz, cp, dcp)
    !
    ! Purpose:
    !
    ! Computes the fourier coefficients of the legendre
    ! polynomial p_n^0 and its derivative.
    ! n is the degree and n/2 or (n+1)/2
    ! coefficients are returned in cp depending on whether
    ! n is even or odd. The same number of coefficients
    ! are returned in dcp. For n even the constant
    ! coefficient is returned in cz.
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)  :: n
    real (wp),    intent (out) :: cz
    real (wp),    intent (out) :: cp(n/2+1)
    real (wp),    intent (out) :: dcp(n/2+1)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer :: j !! Counter
    real    :: t1, t2, t3, t4
    !----------------------------------------------------------------------

    associate( ncp => (n+1)/2 )

        t1 = -1.0_wp
        t2 = real(n + 1, kind=wp)
        t3 = 0.0_wp
        t4 = real(2*n + 1, kind=wp)

        if (mod(n, 2) == 0) then
            cp(ncp) = 1.0_wp
            do j = ncp, 2, -1
                t1 = t1+2.0_wp
                t2 = t2-1.0_wp
                t3 = t3+1.0_wp
                t4 = t4-2.0_wp
                cp(j-1) = (t1*t2)/(t3*t4)*cp(j)
            end do
            t1 = t1+2.0_wp
            t2 = t2-1.0_wp
            t3 = t3+1.0_wp
            t4 = t4-2.0_wp
            cz = (t1*t2)/(t3*t4)*cp(1)
            do j=1, ncp
                dcp(j) = real(2*j, kind=wp)*cp(j)
            end do
        else
            cp(ncp) = 1.0_wp
            do j = ncp-1, 1, -1
                t1 = t1+2.0_wp
                t2 = t2-1.0_wp
                t3 = t3+1.0_wp
                t4 = t4-2.0_wp
                cp(j) = (t1*t2)/(t3*t4)*cp(j+1)
            end do
            do j=1, ncp
                dcp(j) = real(2*j-1, kind=wp)*cp(j)
            end do
        end if
    end associate

end subroutine cpdp


pure subroutine tpdp(n, theta, cz, cp, dcp, pb, dpb)
    !
    ! Purpose:
    !
    ! Computes pn(theta) and its derivative dpb(theta) with
    ! respect to theta
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)  :: n
    real (wp),    intent (in)  :: theta
    real (wp),    intent (in)  :: cz
    real (wp),    intent (in)  :: cp(n/2+1)
    real (wp),    intent (in)  :: dcp(n/2+1)
    real (wp),    intent (out) :: pb
    real (wp),    intent (out) :: dpb
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip) :: k, kdo
    real (wp)    :: chh, cth, sth
    !----------------------------------------------------------------------

    associate( &
        cdt => cos(2.0_wp * theta), &
        sdt => sin(2.0_wp * theta) &
        )

        if (mod(n, 2) == 0) then
            !
            !==> n even
            !
            kdo = n/2
            pb = 0.5_wp * cz
            dpb = 0.0_wp
            if (n > 0) then
                cth = cdt
                sth = sdt
                do k=1, kdo
                    !      pb = pb+cp(k)*cos(2*k*theta)
                    pb = pb+cp(k)*cth
                    !      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
                    dpb = dpb-dcp(k)*sth
                    chh = cdt*cth-sdt*sth
                    sth = sdt*cth+cdt*sth
                    cth = chh
                end do
            end if
        else
            !
            !==> n odd
            !
            kdo = (n + 1)/2
            pb = 0.0_wp
            dpb = 0.0_wp
            cth = cos(theta)
            sth = sin(theta)
            do k=1, kdo
                pb = pb+cp(k)*cth
                dpb = dpb-dcp(k)*sth
                chh = cdt*cth-sdt*sth
                sth = sdt*cth+cdt*sth
                cth = chh
            end do
        end if
    end associate

end subroutine tpdp
