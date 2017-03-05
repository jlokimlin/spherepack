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
!     6/98
!
!     a program for testing all vorticity and ivnerse vorticity routines
!
!
!     (1) first set a stream function and velocity potential scalar fields as
!         polys in x, y, z restricted to the sphere
!
!     (2) derive a vector field (v, w) from (1)
!
!     (3) compute the vorticity vt of (2) and compare with the vorticity
!         computed analytically
!
!     (4) compute vector field (ve, we) using br, bi, cr, ci from (v, w) with
!         br=bi=0.0
!
!     (5) invert the vorticity in (3) and compare with (4)
!

program tvrt
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
    real(wp) :: dwork
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
    integer(ip) :: ldwork
    integer(ip) :: lldwork
    integer(ip) :: lleng
    integer(ip) :: llsav
    integer(ip) :: lsave
    integer(ip) :: lwork
    integer(ip) :: mdab
    integer(ip) :: mdc
    integer(ip) :: mmdab
    integer(ip) :: mmdc
    integer(ip) :: nlat
    integer(ip) :: nlon
    integer(ip) :: nmax
    integer(ip) :: nnlat
    integer(ip) :: nnlon
    integer(ip) :: nnt
    integer(ip) :: nt
    real(wp) :: pertrb
    real(wp) :: phi

    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: thetag
    real(wp) :: v
    real(wp) :: ve
    real(wp) :: vt
    real(wp) :: vte
    real(wp) :: w
    real(wp) :: we
    real(wp) :: work
    real(wp) :: wsave
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 24, nnlon= 14, nnt = 3)
    parameter (mmdab = (nnlon+2)/2, mmdc = (nnlon+1)/2)
    parameter (lleng= 5*nnlat*nnlat*nnlon, llsav=5*nnlat*nnlat*nnlon)
    dimension work(lleng), wsave(llsav)
    parameter (lldwork = 4*nnlat*nnlat)
    dimension dwork(lldwork)
    dimension br(mmdc, nnlat, nnt), bi(mmdc, nnlat, nnt)
    dimension cr(mmdc, nnlat, nnt), ci(mmdc, nnlat, nnt)
    dimension a(mmdab, nnlat, nnt), b(mmdab, nnlat, nnt)
    dimension vt(nnlat, nnlon, nnt)
    dimension thetag(nnlat), dtheta(nnlat), dwts(nnlat)
    dimension v(nnlat, nnlon, nnt), w(nnlat, nnlon, nnt)
    dimension ve(nnlat, nnlon, nnt), we(nnlat, nnlon, nnt)
    dimension pertrb(nnt)
    real dtheta, dwts
    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    nmax = max(nlat, nlon)
    mdab = mmdab
    mdc = mmdc
    lwork = lleng
    lsave = llsav
    ldwork = lldwork
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
    call compute_gaussian_latitudes_and_weights(nlat, dtheta, dwts, ier)
    do  i=1, nlat
        thetag(i) = dtheta(i)
    end do
    call name("gaqd")
    call iout(ier, " ier")
    call vecout(thetag, "thtg", nlat)
    !
    !     test all vorticity subroutines
    !
    do icase=1, 4
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !
        !     set scalar stream and velocity potential fields as polys in x, y, z
        !     and then set v, w from st, sv scalar fields
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
                    dxdp = -sint*sinp
                    dydt = cost*sinp
                    dydp = sint*cosp
                    dzdt = -sint
                    dzdp = 0.0
                    select case (k)
                        case (1)
                            !              st(i, j, k) = x
                            !              sv(i, j, k) = y
                            !
                            !          v = -1/sin(theta)*dstdp + dsvdt
                            !
                            !          w =  1/sin(theta)*dsvdp + dstdt
                            !
                            v(i, j, k) = sinp + cost*sinp
                            w(i, j, k) = cosp + cost*cosp
                            vt(i, j, k) = -2.0*sint*cosp
                        case (2)
                            !              st = y
                            !              sv = z
                            v(i, j, k) = -cosp-sint
                            w(i, j, k) = cost*sinp
                            !         sint*vt = -dvdp + sint*dwdt + cost*w
                            !                 = sinp + sint*(-sint*sinp)+cost*cost*sinp
                            !                 = sinp + (cost**2-sint**2)*sinp
                            vt(i, j, k) = -2.*sint*sinp
                        case (3)
                            !           st = x
                            !           sv = z
                            v(i, j, k) = sinp - sint
                            w(i, j, k) = cost*cosp
                            !     sint*vt = -cosp-sint*sint*sinp+cost*cost*cosp
                            !             = -cosp + (1-2.*sint**2)*cosp =
                            vt(i, j, k) = -2.*sint*cosp
                    end select
                end do
            end do
        end do

        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(v(1, 1, kk), "   v", nlat, nlon)
        !     call aout(w(1, 1, kk), "   w", nlat, nlon)
        !     call aout(vt(1, 1, kk), "  vt", nlat, nlon)
        !     end do

        if (icase==1) then

            call name("**ec")
            !
            !     analyze vector field
            !
            call vhaeci(nlat, nlon, wsave, ierror)
            call name("vhai")
            call iout(ierror, "ierr")

            call vhaec(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdc, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vha ")
            call iout(ierror, "ierr")

            !     if (nmax.lt.10) then
            !     do kk=1, nt
            !     call iout(kk, "**kk")
            !     call aout(br(1, 1, kk), "  br", mdc, nlat)
            !     call aout(bi(1, 1, kk), "  bi", mdc, nlat)
            !     call aout(cr(1, 1, kk), "  cr", mdc, nlat)
            !     call aout(ci(1, 1, kk), "  ci", mdc, nlat)
            !     end do
            !     end if
            !
            !     compute vorticity of (v, w) in vt
            !

            call shseci(nlat, nlon, wsave, ierror)

            call name("vrti")
            call iout(ierror, "ierr")

            call vrtec(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wsave, lsave, work, lwork, ierror)

            call name("vrt ")
            call iout(ierror, "ierr")
            call iout(nlat, "nlat")
            call iout(nlon, "nlon")

        else if (icase==2) then

            call name("**es")
            call shsesi(nlat, nlon, wsave, ierror)

            call name("vrti")
            call iout(ierror, "ierr")

            call vrtes(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wsave, lsave, work, lwork, ierror)

            call name("vrt ")
            call iout(ierror, "ierr")

        else if (icase == 3) then

            call name("**gc")

            call shsgci(nlat, nlon, wsave, ierror)

            call name("vrti")
            call iout(ierror, "ierr")

            call vrtgc(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wsave, lsave, work, lwork, ierror)

            call name("vrt ")
            call iout(ierror, "ierr")

        else if (icase == 4) then

            call name("**gs")

            call shsgsi(nlat, nlon, wsave, ierror)

            call name("vrti")
            call iout(ierror, "ierr")

            call vrtgs(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdc, nlat, &
                wsave, lsave, work, lwork, ierror)

            call name("vrt ")
            call iout(ierror, "ierr")
        end if

        !     if (nmax.lt.10) then
        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(vt(1, 1, kk), "  vt", nlat, nlon)
        !     end do
        !     end if
        !
        !     compute "error" in vt
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
                    select case (k)
                        case (1)
                            vte = -2.0*sint*cosp
                        case (2)
                            vte = -2.*sint*sinp
                        case (3)
                            vte = -2.*sint*cosp
                    end select
                    err2 = err2 + (vt(i, j, k)-vte)**2
                end do
            end do
        end do
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2, "err2")
        !
        !     now recompute (v, w) inverting vt using ivrt(ec, es, gc, gs)
        !     and compare with (ve, we) generated by synthesizing br, bi, cr, ci
        !     with br=bi=0.0
        !

        do k=1, nt
            do i=1, mdc
                do j=1, nlat
                    br(i, j, k) = 0.0
                    bi(i, j, k) = 0.0
                end do
            end do
        end do

        select case (icase)
            case (1)
        		
                call name("**ec")
        		
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call vhseci(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhsec(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wsave, lsave, work, lwork, ierror)
        		
                call name("vhs ")
                call iout(ierror, "ierr")
        		
                call shaeci(nlat, nlon, wsave, ierror)
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shaec(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wsave, ierror)
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
        		
                call vhseci(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call ivrtec(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
                call name("ivrt")
                call iout(ierror, "ierr")
                call vout(pertrb, "prtb")
            case (2)
        		
                call name("**es")
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call vhsesi(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhses(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wsave, lsave, work, lwork, ierror)
                call name("vhs ")
                call iout(ierror, "ierr")
        		
        		
                call shaesi(nlat, nlon, wsave, ierror)
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shaes(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wsave, ierror)
                call name("sha ")
                call iout(ierror, "ierr")
                call iout(lsave, "lsav")
                call iout(lwork, "lwrk")
        		
                call vhsesi(nlat, nlon, wsave, ierror)
                call name("ivti")
                call iout(ierror, "ierr")
        		
                call ivrtes(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
                call name("ivrt")
                call iout(ierror, "ierr")
                call vout(pertrb, "prtb")
            case (3)
        		
                call name("**gc")
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call vhsgci(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhsgc(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wsave, lsave, work, lwork, ierror)
        		
                call name("vhs ")
                call iout(ierror, "ierr")
        		
                call shagci(nlat, nlon, wsave, ierror)
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shagc(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wsave, ierror)
                call name("sha ")
                call iout(ierror, "ierr")
                call iout(lsave, "lsav")
                call iout(lwork, "lwrk")
        		
                call vhsgci(nlat, nlon, wsave, ierror)
                call name("ivti")
                call iout(ierror, "ierr")
        		
                call ivrtgc(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
        		
                call name("ivrt")
                call iout(ierror, "ierr")
                call vout(pertrb, "prtb")
            case (4)
        		
                call name("**gs")
                !
                !     set vector field (ve, we) with br=bi=0.0 for comparison with inverted vt
                !
                call vhsgsi(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhsgs(nlat, nlon, ityp, nt, ve, we, nlat, nlon, br, bi, cr, ci, &
                    mdc, nlat, wsave, lsave, work, lwork, ierror)
                call name("vhs ")
                call iout(ierror, "ierr")
        		
                call shagsi(nlat, nlon, wsave, ierror)
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shags(nlat, nlon, isym, nt, vt, nlat, nlon, a, b, &
                    mdab, nlat, wsave, ierror)
                call name("sha ")
                call iout(ierror, "ierr")
                call iout(lsave, "lsav")
                call iout(lwork, "lwrk")
        		
                call vhsgsi(nlat, nlon, wsave, ierror)
                call name("ivti")
                call iout(ierror, "ierr")
        		
                call ivrtgs(nlat, nlon, isym, nt, v, w, nlat, nlon, a, b, &
                    mdab, nlat, wsave, lsave, work, lwork, pertrb, ierror)
                call name("ivrt")
                call iout(ierror, "ierr")
                call vout(pertrb, "prtb")
        end select


        !     if (nmax.lt.10) then
        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(v(1, 1, kk), "   v", nlat, nlon)
        !     call aout(w(1, 1, kk), "   w", nlat, nlon)
        !     end do
        !     end if

        !
        !     compare this v, w with ve, we
        !
        err2v = 0.0
        err2w = 0.0
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    err2v = err2v + (v(i, j, k)-ve(i, j, k))**2
                    err2w = err2w + (w(i, j, k)-we(i, j, k))**2
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

end program tvrt
