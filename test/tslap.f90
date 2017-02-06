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
!     3/97
!
!     a program for testing slap, islap (ec, es, gc, gs)
!
!     (1) set a scalar field s as poly in x, y, z restricted to sphere
!
!     (2) compute scalar laplacian in array sclp using slap(ec, es, gc, gs)
!
!     (3) compare (2) with analytic scalar laplacian in sclpe
!
!     (4) compute the inverse  of (2) and compare with (1)
!
program tslap
use spherepack
    implicit none
    real(wp) :: a
    real(wp) :: b
    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: d2xdp2
    real(wp) :: d2xdt2
    real(wp) :: d2ydp2
    real(wp) :: d2ydt2
    real(wp) :: d2zdp2
    real(wp) :: d2zdt2
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: dxdp
    real(wp) :: dxdt
    real(wp) :: dydp
    real(wp) :: dydt
    real(wp) :: dzdp
    real(wp) :: dzdt
    real(wp) :: err2
    real(wp) :: err2s
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: ier
    integer(ip) :: ierror
    integer(ip) :: isym
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: ldwork
    integer(ip) :: lldwork
    integer(ip) :: lleng
    integer(ip) :: llsav
    integer(ip) :: lsave
    integer(ip) :: lwork
    integer(ip) :: mdab
    integer(ip) :: mmdab
    integer(ip) :: nlat
    integer(ip) :: nlon
    integer(ip) :: nnlat
    integer(ip) :: nnlon
    integer(ip) :: nnt
    integer(ip) :: nt
    real(wp) :: phi

    real(wp) :: ptrb
    real(wp) :: s
    real(wp) :: sclp
    real(wp) :: sclpe
    real(wp) :: se
    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: thetag
    real(wp) :: work
    real(wp) :: wsave
    real(wp) :: x
    real(wp) :: xlm
    real(wp) :: y
    real(wp) :: z
    parameter(nnlat=15, nnlon= 22, nnt = 3)
    parameter (mmdab = (nnlon+2)/2)

    parameter (lleng=15*nnlat*nnlat*nnlon, llsav=5*nnlat*nnlat*nnlon)
    dimension work(lleng), wsave(llsav)
    parameter (lldwork = nnlat*(nnlat+4))
    real dwork(lldwork)
    dimension a(mmdab, nnlat, nnt), b(mmdab, nnlat, nnt), s(nnlat, nnlon, nnt)
    dimension sclp(nnlat, nnlon, nnt)
    dimension sclpe(nnlat, nnlon, nnt)
    dimension ptrb(nnt), xlm(nnt)
    dimension thetag(nnlat), dtheta(nnlat), dwts(nnlat)
    real dtheta, dwts
    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    mdab = mmdab
    lwork = lleng
    lsave = llsav
    nt = nnt
    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")
    isym = 0
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
    !     set helmholtz constant zero for laplacian inversion
    !
    do k=1, nt
        xlm(k) = 0.0
    end do
    !
    !     test all analysis and synthesis subroutines
    !
    do icase=1, 4
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !
        !     set scalar field as poly in x, y, z restricted to the sphere
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
                    d2xdt2 = -x
                    dxdp = -cost*sinp
                    d2xdp2 = -x
                    dydt = cost*sinp
                    d2ydt2 = -y
                    dydp = cost*cosp
                    d2ydp2 = -y
                    dzdt = -sint
                    d2zdt2 = -z
                    dzdp = 0.
                    d2zdp2 = 0.
                    select case (k)
                    	case (1)
                    		s(i, j, k) = x+y
                    		sclpe(i, j, k) = -2.*(x+y)
                    	case (2)
                    		s(i, j, k) = x+z
                    		sclpe(i, j, k) = -2.*(x+z)
                    	case (3)
                    		s(i, j, k) = y+z
                    		sclpe(i, j, k) = -2.*(y+z)
                    end select
                end do
            end do
        end do
        !     do k=1, nt
        !     call iout(k, "   k")
        !     call aout(s(1, 1, k), "   s", nlat, nlon)
        !     call aout(sclpe(1, 1, k), "sclp", nlat, nlon)
        !     end do

        !     call aout(s, "   s", nlat, nlon)


        select case (icase)
        	case (1)
        		
        		call name("**ec")
        		
        		call shaeci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shaec(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wsave, &
        		lsave, work, lwork, ierror)
        		call name("sha ")
        		call iout(ierror, "ierr")
        		
        		call shseci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call slapec(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		call name("slap")
        		call iout(ierror, "ierr")
        	case (2)
        		
        		call name("**es")
        		
        		call shaesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shaes(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wsave, &
        		lsave, work, lwork, ierror)
        		call name("sha ")
        		call iout(ierror, "ierr")
        		
        		call shsesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call slapes(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		call name("slap")
        		call iout(ierror, "ierr")
        	case (3)
        		
        		call name("**gc")
        		
        		call shagci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shagc(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wsave, &
        		lsave, work, lwork, ierror)
        		call name("sha ")
        		call iout(ierror, "ierr")
        		
        		call shsgci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call slapgc(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		call name("slap")
        		call iout(ierror, "ierr")
        	case (4)
        		
        		call name("**gs")
        		
        		call shagsi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call name("shai")
        		call iout(ierror, "ierr")
        		
        		call shags(nlat, nlon, isym, nt, s, nlat, nlon, a, b, mdab, nlat, wsave, &
        		lsave, work, lwork, ierror)
        		call name("sha ")
        		call iout(ierror, "ierr")
        		
        		call shsgsi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call slapgs(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		call name("slap")
        		call iout(ierror, "ierr")
        end select
        !
        !     compute "error" in sclp
        !
        err2 = 0.0
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    err2 = err2 + (sclpe(i, j, k)-sclp(i, j, k))**2
                end do
            end do
        !     call iout(k, "   k")
        !     call aout(sclp(1, 1, k), "sclp", nlat, nlon)
        end do
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2, "err2")
        !
        !     invert sclp
        !
        select case (icase)
        	case (1)
        		
        		call shaeci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call shaec(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		
        		call shseci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call islapec(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ptrb, ierror)
        		call name("isla")
        		call iout(ierror, "ierr")
        		call vecout(ptrb, "ptrb", nt)
        	case (2)
        		
        		call shaesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call shaes(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		
        		call shsesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call islapes(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ptrb, ierror)
        		call name("isla")
        		call iout(ierror, "ierr")
        		call vecout(ptrb, "ptrb", nt)
        	case (3)
        		
        		call shagci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call shagc(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		
        		call shsgci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call islapgc(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ptrb, ierror)
        		call name("isla")
        		call iout(ierror, "ierr")
        		call vecout(ptrb, "ptrb", nt)
        	case (4)
        		
        		call name("**gs")
        		
        		call shagsi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call shags(nlat, nlon, isym, nt, sclp, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ierror)
        		
        		call shsgsi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		call name("shsi")
        		call iout(ierror, "ierr")
        		
        		call islapgs(nlat, nlon, isym, nt, xlm, s, nlat, nlon, a, b, mdab, nlat, &
        		wsave, lsave, work, lwork, ptrb, ierror)
        		call name("isla")
        		call iout(ierror, "ierr")
        		call vecout(ptrb, "ptrb", nt)
        end select

        !     call aout(s, "   s", nlat, nlon)


        !
        !     compare s with original
        !
        err2s = 0.0
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
                    select case (k)
                    	case (1)
                    		se = x+y
                    	case (2)
                    		se = x+z
                    	case (3)
                    		se = y+z
                    end select
                    err2s = err2s+(s(i, j, k) - se)**2
                end do
            end do
        end do
        err2s = sqrt(err2s/(nlat*nlon*nt))
        call vout(err2s, "errs")
    !
    !     end of icase loop
    !
    end do

end program tslap
