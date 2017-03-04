!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK                               *
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
!
!     12/96
!
!     a program for testing all divergence and inverse divergence routines
!
!     (1) first set an irrotational vector field v, w
!
!     (2) compute the coefficients br, bi, cr, ci of (v, w) using vector analysis
!
!     (3) compute the divergence of (v, w) using divec, dives, divgc, divgs
!
!     (4) analystically compute the divergence and compare with (3)
!
!     (5) invert the divergence and compare with the irrotational (v, w)
!
program tdiv
use spherepack
    implicit none
    real(wp) :: a
    real(wp) :: b
    real(wp) :: bi
    real(wp) :: br
    real(wp) :: ci
    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: cr
    real(wp) :: d2xdp2
    real(wp) :: d2xdt2
    real(wp) :: d2ydp2
    real(wp) :: d2ydt2
    real(wp) :: d2zdp2
    real(wp) :: d2zdt2
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: dsfdp
    real(wp) :: dsfdt
    real(wp) :: dv
    real(wp) :: dvdt
    real(wp) :: dve
    real(wp) :: dwork
    real(wp) :: dx2dt2
    real(wp) :: dxdp
    real(wp) :: dxdt
    real(wp) :: dydp
    real(wp) :: dydt
    real(wp) :: dzdp
    real(wp) :: dzdt
    real(wp) :: err2
    real(wp) :: err2v
    real(wp) :: err2w
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: ier
    integer(ip) :: ierror
    integer(ip) :: isym
    integer(ip) :: ityp
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: kk
    integer(ip) :: ldwork
    integer(ip) :: lldwork
    integer(ip) :: lleng
    integer(ip) :: llsav
    integer(ip) :: lsave
    integer(ip) :: lwork
    integer(ip) :: mdab
    integer(ip) :: mdb
    integer(ip) :: mmdab
    integer(ip) :: mmdb
    integer(ip) :: nlat
    integer(ip) :: nlon
    integer(ip) :: nmax
    integer(ip) :: nnlat
    integer(ip) :: nnlon
    integer(ip) :: nnt
    integer(ip) :: nt
    real(wp) :: pertrb
    real(wp) :: phi

    real(wp) :: sf
    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: thetag
    real(wp) :: v
    real(wp) :: ve
    real(wp) :: w
    real(wp) :: we
    real(wp) :: work
    real(wp) :: wsave
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    real(wp) :: dtheta, dwts
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 15, nnlon=  9, nnt = 2)
    ! *** be sure to set
    !     mmdb = min(nlat, (nlon+1)/2), mmdab = min(nlat, (nnlon+2)/2)
    parameter (mmdab = (nnlon+2)/2, mmdb = (nnlon+1)/2)
    parameter (lleng= 5*nnlat*nnlat*nnlon, llsav=15*nnlat*nnlat*nnlon)
    parameter (lldwork = 4*nnlat*nnlat)
    dimension work(lleng), wsave(llsav)
    dimension br(mmdb, nnlat, nnt), bi(mmdb, nnlat, nnt)
    dimension dwork(lldwork)
    dimension cr(mmdb, nnlat, nnt), ci(mmdb, nnlat, nnt)
    dimension a(mmdab, nnlat, nnt), b(mmdab, nnlat, nnt)
    dimension dv(nnlat, nnlon, nnt)
    dimension thetag(nnlat), dtheta(nnlat), dwts(nnlat)
    dimension v(nnlat, nnlon, nnt), w(nnlat, nnlon, nnt)
    dimension pertrb(nnt)

    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    nmax = max(nlat, nlon)
    mdab = mmdab
    mdb = mmdb

    lwork = lleng
    lsave = llsav
    nt = nnt
    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")
    isym = 0
    ityp = 0
    !
    !     set equally spaced colatitude and longitude increments
    !
    dphi = (pi+pi)/nlon
    dlat = pi/(nlat-1)
    !
    !     compute nlat gaussian points in thetag
    !
    ldwork = lldwork
    call compute_gaussian_latitudes_and_weights(nlat, dtheta, dwts, ier)
    do  i=1, nlat
        thetag(i) = dtheta(i)
    end do
    call name("gaqd")
    call iout(ier, " ier")
    call vecout(thetag, "thtg", nlat)
    !
    !     test all divergence and inverse divergence subroutines
    !
    do icase=1, 4
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !     set vector field v, w
        !
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=thetag(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    d2xdt2 = -sint*cosp
                    dxdp = -sint*sinp
                    d2xdp2 = -sint*cosp
                    dydt = cost*sinp
                    d2ydt2 = -sint*sinp
                    dydp = sint*cosp
                    d2ydp2 = -sint*sinp
                    dzdt = -sint
                    d2zdt2 = -cost
                    dzdp = 0.0
                    d2zdp2 = 0.0
                    if (k==1) then
                        sf = x*y
                        dsfdt = x*dydt+y*dxdt
                        dsfdp = x*dydp+y*dxdp
                        v(i, j, k) = dsfdt
                        w(i, j, k) = cosp*dydp+sinp*dxdp
                        !              dv = 1/sint*(d(sint*v)/dt + dw/dp)
                        dvdt = x*d2ydt2 + 2.*dxdt*dydt + y*d2xdt2
                        !              1/sint*dwdp = 1/sint*(cosp*d2ydp2-sinp*dydp+sinp*d2xdp2+cosp*dxdp)
                        !                          = -4.*sinp*cosp
                        dv(i, j, k) = dvdt + cost*(cosp*dydt+sinp*dxdt) -4.*cosp*sinp
                    else if (k==2) then
                        sf = x*z
                        dsfdt = x*dzdt+z*dxdt
                        dsfdp = x*dzdp+z*dxdp
                        v(i, j, k) = dsfdt
                        w(i, j, k) = -cost*sinp
                        dvdt = x*d2zdt2+2.*dzdt*dxdt + z*dx2dt2
                        !              dv = 1/sint*(d(sint*v)/dt + dw/dp)
                        dv(i, j, k) = -6.*cost*sint*cosp
                    else if (k==3) then
                        sf = y*z
                        dsfdt = y*dzdt+z*dydt
                        dsfdp = y*dzdp+z*dydp
                        v(i, j, k) = dsfdt
                        w(i, j, k) = z*cosp
                        dv(i, j, k) = -6.*cost*sint*sinp
                    end if
                end do
            end do
        end do

        !     if (nmax.lt.10) then
        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(v(1, 1, kk), "   v", nlat, nlon)
        !     call aout(w(1, 1, kk), "   w", nlat, nlon)
        !     call aout(dv(1, 1, kk), "  dv", nlat, nlon)
        !     end do
        !     end if

        if (icase==1) then

            call name("**ec")
            !
            !     analyze vector field
            !
            call vhaeci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
            call name("shs ")
            call iout(ierror, "ierr")

            call vhaec(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdb, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vha ")
            call iout(ierror, "ierr")

            !c    if (nmax.lt.10) then
            !c    do kk=1, nt
            !c    call iout(kk, "**kk")
            !c    call aout(br(1, 1, kk), "  br", mdb, nlat)
            !c    call aout(bi(1, 1, kk), "  bi", mdb, nlat)
            !c    call aout(cr(1, 1, kk), "  cr", mdb, nlat)
            !c    call aout(ci(1, 1, kk), "  ci", mdb, nlat)
            !c    end do
            !c    end if
            !
            !     compute divergence of (v, w) in dv
            !

            call shseci(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call divec(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("div ")
            call iout(ierror, "ierr")
            call iout(nlat, "nlat")
            call iout(nlon, "nlon")


        else if (icase==2) then

            call name("**es")
            call shsesi(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call dives(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)

            call name("div ")
            call iout(ierror, "ierr")

        else if (icase == 3) then

            call name("**gc")

            call shsgci(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call divgc(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("div ")
            call iout(ierror, "ierr")

        else if (icase == 4) then

            call name("**gs")

            call shsgsi(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call divgs(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("div ")
            call iout(ierror, "ierr")

        end if

        !     if (nmax.lt.10) then
        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(dv(1, 1, kk), "  dv", nlat, nlon)
        !     end do
        !     end if
        !
        !     compute "error" in dv
        !
        err2 = 0.0
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=thetag(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    d2xdt2 = -sint*cosp
                    dxdp = -sint*sinp
                    d2xdp2 = -sint*cosp
                    dydt = cost*sinp
                    d2ydt2 = -sint*sinp
                    dydp = sint*cosp
                    d2ydp2 = -sint*sinp
                    dzdt = -sint
                    d2zdt2 = -cost
                    dzdp = 0.0
                    d2zdp2 = 0.0
                    if (k==1) then
                        sf = x*y
                        dvdt = x*d2ydt2 + 2.*dxdt*dydt + y*d2xdt2
                        dve = dvdt + cost*(cosp*dydt+sinp*dxdt) - 4.*cosp*sinp
                    else if (k==2) then
                        !              sf = x*z
                        dve = -6.*sint*cost*cosp
                    else if (k==3) then
                        !              sf = y*z
                        dve = -6.*cost*sint*sinp
                    end if
                    err2 = err2 + (dv(i, j, k)-dve)**2
                end do
            end do
        end do
        !
        !     set and print least squares error in dv
        !
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2, "err2")
        !
        !     now recompute (v, w) inverting dv using idiv(ec, es, gc, gs)
        !
        do kk=1, nt
            do j=1, nlon
                do i=1, nlat
                    v(i, j, kk) = 0.0
                    w(i, j, kk) = 0.0
                end do
            end do
        end do


        if (icase==1) then

            call name("**ec")


            call shaeci(nlat, nlon, wsave, ierror)
            call name("shai")
            call iout(ierror, "ierr")

            call shaec(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, ierror)
            call name("sha ")
            call iout(ierror, "ierr")
            call iout(lsave, "lsav")
            call iout(lwork, "lwrk")

            !     if (nmax.lt.10) then
            !     do kk=1, nt
            !     call iout(kk, "**kk")
            !     call aout(a(1, 1, kk), "   a", nlat, nlat)
            !     call aout(b(1, 1, kk), "   b", nlat, nlat)
            !     end do
            !     end if

            call vhseci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
            call name("idvi")
            call iout(ierror, "ierr")

            call idivec(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
            call name("idiv")
            call iout(ierror, "ierr")
            call vecout(pertrb, "prtb", nt)

        else if (icase==2) then

            call name("**es")


            call shaesi(nlat, nlon, wsave, ierror)
            call name("shai")
            call iout(ierror, "ierr")

            call shaes(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, ierror)
            call name("sha ")
            call iout(ierror, "ierr")
            call iout(lsave, "lsav")
            call iout(lwork, "lwrk")

            call vhsesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
            call name("idvi")
            call iout(ierror, "ierr")

            call idives(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
            call name("idiv")
            call iout(ierror, "ierr")
            call vecout(pertrb, "prtb", nt)

        else if (icase==3) then

            call name("**gc")


            call shagci(nlat, nlon, wsave, ierror)
            call name("shai")
            call iout(ierror, "ierr")

            call shagc(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, ierror)
            call name("sha ")
            call iout(ierror, "ierr")
            call iout(lsave, "lsav")
            call iout(lwork, "lwrk")

            call vhsgci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
            call name("idvi")
            call iout(ierror, "ierr")

            call idivgc(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
            call name("idiv")
            call iout(ierror, "ierr")
            call vecout(pertrb, "prtb", nt)

        else if (icase==4) then

            call name("**gs")


            call shagsi(nlat, nlon, wsave, ierror)
            call name("shai")
            call iout(ierror, "ierr")

            call shags(nlat, nlon, isym, nt, dv, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, ierror)
            call name("sha ")
            call iout(ierror, "ierr")
            call iout(lsave, "lsav")
            call iout(lwork, "lwrk")

            call vhsgsi(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
            call name("idvi")
            call iout(ierror, "ierr")

            call idivgs(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
            call name("idiv")
            call iout(ierror, "ierr")
            call vecout(pertrb, "prtb", nt)

        end if


        !     if (nmax.lt.10) then
        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(v(1, 1, kk), "   v", nlat, nlon)
        !     call aout(w(1, 1, kk), "   w", nlat, nlon)
        !     end do
        !     end if

        !
        !     compare this v, w with original
        !
        err2v = 0.0
        err2w = 0.0
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=thetag(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    dxdp = -sint*sinp
                    dydt = cost*sinp
                    dydp = sint*cosp
                    dzdt = -sint
                    dzdp = 0.0
                    if (k==1) then
                        sf = x*y
                        dsfdt = x*dydt+y*dxdt
                        dsfdp = x*dydp+y*dxdp
                        ve = dsfdt
                        we = cosp*dydp+sinp*dxdp
                    else if (k==2) then
                        sf = x*z
                        dsfdt = x*dzdt+z*dxdt
                        dsfdp = x*dzdp+z*dxdp
                        ve = dsfdt
                        we = -cost*sinp
                    else if (k==3) then
                        sf = y*z
                        dsfdt = y*dzdt+z*dydt
                        dsfdp = y*dzdp+z*dydp
                        ve = dsfdt
                        we = z*cosp
                    else if (k==4) then
                    end if
                    err2v = err2v + (v(i, j, k)-ve)**2
                    err2w = err2w + (w(i, j, k)-we)**2
                end do
            end do
        end do
        err2v = sqrt(err2v/(nlat*nlon*nt))
        err2w = sqrt(err2w/(nlat*nlon*nt))
        call vout(err2v, "errv")
        call vout(err2w, "errw")

    !
    !     end of icase loop
    !
    end do

end program tdiv

