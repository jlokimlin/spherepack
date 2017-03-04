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
!     3/97
!
!     a program for testing all divergence, vorticity and idvt(ec, es, gc, gs) routines
!
!     (1) first set a valid vector field by setting a stream function sf and velocity
!         potential function sv as polys in x, y, z restricted to the sphere.  Then
!         derive (v, w) and dv, vt from sf and sv analytically by differentiation.
!         (see tvha.f)
!
!     (2) compute the coefficients br, bi, cr, ci of (v, w) using vector analysis
!
!     (3) compute the divergence and vorticity of (v, w) using div, vrt (es, ec, gc, gs)
!
!     (4) compare with divergence and vorticity from (1)
!
!     (5) invert dv, vt with idvt(ec, es, gc, gs) and compare with vector field from (1)
!
program tidvt
use spherepack
    implicit none
    real(wp) :: ad
    real(wp) :: av
    real(wp) :: bd
    real(wp) :: bi
    real(wp) :: br
    real(wp) :: bv
    real(wp) :: ci
    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: cr
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: dv
    real(wp) :: dvsav
    real(wp) :: err2d
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
    real(wp) :: phi

    real(wp) :: prtbd
    real(wp) :: prtbv
    real(wp) :: ptrbd
    real(wp) :: ptrbv
    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: thetag
    real(wp) :: v
    real(wp) :: vsav
    real(wp) :: vt
    real(wp) :: vtsav
    real(wp) :: w
    real(wp) :: work
    real(wp) :: wsav
    real(wp) :: wsave
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 25, nnlon= 16, nnt = 3)
    parameter(mmdab = (nnlon+2)/2, mmdb = (nnlon+1)/2)
    parameter (lleng= 5*nnlat*nnlat*nnlon, llsav=15*nnlat*nnlat*nnlon)
    dimension work(lleng), wsave(llsav)
    parameter (lldwork = 4*nnlat*nnlat )
    real dwork(lldwork)
    dimension br(mmdb, nnlat, nnt), bi(mmdb, nnlat, nnt)
    dimension cr(mmdb, nnlat, nnt), ci(mmdb, nnlat, nnt)
    dimension ad(mmdab, nnlat, nnt), bd(mmdab, nnlat, nnt)
    dimension av(mmdab, nnlat, nnt), bv(mmdab, nnlat, nnt)
    dimension dv(nnlat, nnlon, nnt), vt(nnlat, nnlon, nnt)
    dimension dvsav(nnlat, nnlon, nnt), vtsav(nnlat, nnlon, nnt)
    dimension thetag(nnlat), dtheta(nnlat), dwts(nnlat)
    dimension v(nnlat, nnlon, nnt), w(nnlat, nnlon, nnt)
    dimension vsav(nnlat, nnlon, nnt), wsav(nnlat, nnlon, nnt)
    dimension ptrbd(nnt), ptrbv(nnt)
    real dtheta, dwts
    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    mdab = mmdab
    mdb = mmdb
    nmax = max(nlat, nlon)

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
        !
        !     icase=1 corresponds to "ec"
        !     icase=2 corresponds to "es"
        !     icase=3 corresponds to "gc"
        !     icase=4 corresponds to "gs"

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
                    if (k ==1) then
                        !              sf = x
                        !              sv = y
                        !
                        !          v = -1/sint*dstdp + dsvdt
                        !
                        !          w =  1/sint*dsvdp + dstdt
                        !
                        !          dv = 1/sint*[d(sint*v)/dt + dwdp]  = dvdt + ct/st*v + 1/st*dwdp
                        !
                        !          vt = 1/sint*[-dv/dp + d(sint*w)/dt) = dwdt + ct/st*w - 1/st*dvdp

                        v(i, j, k) = sinp + cost*sinp
                        w(i, j, k) = cosp + cost*cosp
                        dv(i, j, k) = -2.0*sint*sinp
                        vt(i, j, k) = -2.0*sint*cosp
                    else if (k==2) then
                        !              sf = y
                        !              sv = z
                        v(i, j, k) = -cosp-sint
                        w(i, j, k) = cost*sinp
                        vt(i, j, k) = -2.*sint*sinp
                        dv(i, j, k) = -2.*cost
                    else if (k==3) then
                        !              st = x
                        !              sv = z
                        v(i, j, k) = sinp - sint
                        w(i, j, k) = cost*cosp
                        vt(i, j, k) = -2.*sint*cosp
                        dv(i, j, k) = -2.*cost
                    end if
                    !
                    !      save derived vector field, vorticity, divergence for latter comparison
                    !
                    vtsav(i, j, k) = vt(i, j, k)
                    dvsav(i, j, k) = dv(i, j, k)
                    vsav(i, j, k) = v(i, j, k)
                    wsav(i, j, k) = w(i, j, k)
                end do
            end do
        end do

        if (icase==1) then

            call name("**ec")
            !
            !     analyze vector field
            !
            call vhaeci(nlat, nlon, wsave, lsave, work, lwork, ierror)
            call name("vhai")
            call iout(ierror, "ierr")

            call vhaec(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdb, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vha ")
            call iout(ierror, "ierr")
            !
            !     compute divergence, vorticity of (v, w) in dv, vt
            !

            call shseci(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call divec(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("div ")
            call iout(ierror, "ierr")

            call vrtec(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("vrt ")
            call iout(ierror, "ierr")


        else if (icase==2) then

            call name("**es")
            !
            !     analyze vector field
            !
            call vhaesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
            call name("vhai")
            call iout(ierror, "ierr")
            call vhaes(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdb, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vha ")
            call iout(ierror, "ierr")

            call shsesi(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call dives(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("div ")
            call iout(ierror, "ierr")
            call vrtes(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("vrt ")
            call iout(ierror, "ierr")

        else if (icase == 3) then

            call name("**gc")
            !
            !     analyze vector field
            !
            call vhagci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
            call name("vhai")
            call iout(ierror, "ierr")
            call vhagc(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdb, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vha ")
            call iout(ierror, "ierr")

            call shsgci(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call divgc(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("div ")
            call iout(ierror, "ierr")
            call vrtgc(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("vrt ")
            call iout(ierror, "ierr")

        else if (icase == 4) then

            call name("**gs")
            !
            !     analyze vector field
            !
            call vhagsi(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
            call name("vhai")
            call iout(ierror, "ierr")
            call vhags(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdb, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vha ")
            call iout(ierror, "ierr")

            call shsgsi(nlat, nlon, wsave, ierror)
            call name("shsi")
            call iout(ierror, "ierr")

            call divgs(nlat, nlon, isym, nt, dv, nlat, nlon, br, bi, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("div ")
            call iout(ierror, "ierr")
            call vrtgs(nlat, nlon, isym, nt, vt, nlat, nlon, cr, ci, mdb, nlat, &
                wsave, lsave, work, lwork, ierror)
            call name("vrt ")
            call iout(ierror, "ierr")

        end if
        !
        !     compute "error" in dv, vt
        !
        err2d = 0.0
        err2v = 0.0
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    err2d = err2d + (dv(i, j, k)-dvsav(i, j, k))**2
                    err2v = err2v + (vt(i, j, k)-vtsav(i, j, k))**2
                end do
            end do
        end do
        !
        !     set and print least squares error in dv
        !
        err2d = sqrt(err2d/(nt*nlat*nlon))
        call vout(err2d, "errd")
        err2v = sqrt(err2v/(nt*nlat*nlon))
        call vout(err2v, "errv")
        !
        !     now recompute (v, w) inverting dv, vt using idvt(ec, es, gc, gs)
        !
        do kk=1, nt
            do j=1, nlon
                do i=1, nlat
                    v(i, j, kk) = 0.0
                    w(i, j, kk) = 0.0
                end do
            end do
        end do


        select case (icase)
        	case (1)
        		
        		call name("**ec")
        		
        		
        		call shaeci(nlat, nlon, wsave, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shaec(nlat, nlon, isym, nt, dv, nlat, nlon, ad, bd, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		call shaec(nlat, nlon, isym, nt, vt, nlat, nlon, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		call name("sha ")
        		call iout(ierror, "ierr")
        		call iout(lsave, "lsav")
        		call iout(lwork, "lwrk")
        		
        		call vhseci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("idvi")
        		call iout(ierror, "ierr")
        		
        		call idvtec(nlat, nlon, isym, nt, v, w, nlat, nlon, ad, bd, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ptrbd, ptrbv, ierror)
        		call name("idvt")
        		call iout(ierror, "ierr")
        		call vecout(prtbd, "prtd", nt)
        		call vecout(prtbv, "prtv", nt)
        	case (2)
        		
        		call name("**es")
        		
        		call shaesi(nlat, nlon, wsave, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shaes(nlat, nlon, isym, nt, dv, nlat, nlon, ad, bd, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		call shaes(nlat, nlon, isym, nt, vt, nlat, nlon, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		
        		call name("sha ")
        		call iout(ierror, "ierr")
        		call iout(lsave, "lsav")
        		call iout(lwork, "lwrk")
        		
        		call vhsesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call name("idvi")
        		call iout(ierror, "ierr")
        		
        		call idvtes(nlat, nlon, isym, nt, v, w, nlat, nlon, ad, bd, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ptrbd, ptrbv, ierror)
        		call name("idvt")
        		call iout(ierror, "ierr")
        		call vecout(prtbd, "prtd", nt)
        		call vecout(prtbv, "prtv", nt)
        	case (3)
        		
        		call name("**gc")
        		
        		call shagci(nlat, nlon, wsave, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shagc(nlat, nlon, isym, nt, dv, nlat, nlon, ad, bd, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		call shagc(nlat, nlon, isym, nt, vt, nlat, nlon, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		
        		call name("sha ")
        		call iout(ierror, "ierr")
        		call iout(lsave, "lsav")
        		call iout(lwork, "lwrk")
        		
        		call vhsgci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("idvi")
        		call iout(ierror, "ierr")
        		
        		call idvtgc(nlat, nlon, isym, nt, v, w, nlat, nlon, ad, bd, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ptrbd, ptrbv, ierror)
        		call name("idvt")
        		call iout(ierror, "ierr")
        		call vecout(prtbd, "prtd", nt)
        		call vecout(prtbv, "prtv", nt)
        	case (4)
        		
        		call name("**gs")
        		
        		call shagsi(nlat, nlon, wsave, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shags(nlat, nlon, isym, nt, dv, nlat, nlon, ad, bd, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		call shags(nlat, nlon, isym, nt, vt, nlat, nlon, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ierror)
        		
        		call name("sha ")
        		call iout(ierror, "ierr")
        		call iout(lsave, "lsav")
        		call iout(lwork, "lwrk")
        		
        		call vhsgsi(nlat, nlon, wsave, ierror)
        		call name("idvi")
        		call iout(ierror, "ierr")
        		
        		call idvtgs(nlat, nlon, isym, nt, v, w, nlat, nlon, ad, bd, av, bv, &
        		mdab, nlat, wsave, lsave, work, lwork, ptrbd, ptrbv, ierror)
        		call name("idvt")
        		call iout(ierror, "ierr")
        		call vecout(prtbd, "prtd", nt)
        		call vecout(prtbv, "prtv", nt)
        end select

        !
        !     compare this v, w with original derived from sf, sv
        !
        err2v = 0.0
        err2w = 0.0
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    err2v = err2v + (v(i, j, k)-vsav(i, j, k))**2
                    err2w = err2w + (w(i, j, k)-wsav(i, j, k))**2
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

end program tidvt

