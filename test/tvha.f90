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
!     11/96
!
!     a program for testing all vector analysis and synthesis subroutines
!
!     (1) first a scalar stream function and a scalar velocity potential function
!         are set in st, sv by restricting polys in x, y, z to the sphere surface
!
!     (2) the vector vield (v, w) is set by analytically differenting the scalar fields in (1)
!         using the standard formula relating a vector field and the stream and velocity
!         potential scalar fields in colatitude X longitude spherical coordinates
!
!          v = -1/sin(theta)*d(st)/dphi + d(sv)/dtheta
!
!          w =  1/sin(theta)*d(sv)/dphi + d(st)/dtheta
!
!     (3) a vector analysis is performed on (v, w)
!
!     (4) a vector synthesis is performed using coeffs from (3)
!
!     (5) the synthesized vector field from (4) is compared with the vector field from (2)
!
!     note:  vhaec, vhaes, vhagc, vhags, vhsec, vhses, vhsgc, vhsgs are all tested!
!
program tvha
    use spherepack
    implicit none
    real(wp) :: bi
    real(wp) :: br
    real(wp) :: ci
    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: cr
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: dstdp
    real(wp) :: dstdt
    real(wp) :: dsvdp
    real(wp) :: dsvdt
    real(wp) :: dxdp
    real(wp) :: dxdt
    real(wp) :: dydp
    real(wp) :: dydt
    real(wp) :: dzdp
    real(wp) :: dzdt
    real(wp) :: err2v
    real(wp) :: err2w
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: ier
    integer(ip) :: ierror
    integer(ip) :: ityp
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: ldwork
    integer(ip) :: lldwork
    integer(ip) :: lleng
    integer(ip) :: llsav
    integer(ip) :: lsave
    integer(ip) :: lwork
    integer(ip) :: nlat
    integer(ip) :: nlon
    integer(ip) :: nnlat
    integer(ip) :: nnlon
    integer(ip) :: nnt
    integer(ip) :: nt
    real(wp) :: phi

    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: st
    real(wp) :: sv
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
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 25, nnlon= 19, nnt = 2)
    parameter (lleng= 5*nnlat*nnlat*nnlon, llsav= 5*nnlat*nnlat*nnlon)
    parameter (lldwork = 4*nnlat*nnlat)
    real dwork(lldwork)
    dimension work(lleng), wsave(llsav)
    dimension br(nnlat, nnlat, nnt), bi(nnlat, nnlat, nnt)
    dimension cr(nnlat, nnlat, nnt), ci(nnlat, nnlat, nnt)
    dimension st(nnlat, nnlon, nnt), sv(nnlat, nnlon, nnt)
    dimension thetag(nnlat), dtheta(nnlat), dwts(nnlat)
    dimension v(nnlat, nnlon, nnt), w(nnlat, nnlon, nnt)
    real dtheta, dwts
    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    lwork = lleng
    lsave = llsav
    nt = nnt
    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")
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
    !     test all analysis and synthesis subroutines
    !
    do icase=1, 4
        !
        !     icase=1 test vhaec, vhsec
        !     icase=2 test vhaes, vhses
        !     icase=3 test vhagc, vhsgc
        !     icase=4 test vhags, vhsgs
        !
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
                    if (icase==3 .or. icase==4) theta = thetag(i)
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
                            st(i, j, k) = x*y
                            sv(i, j, k) = y*z
                            dstdt = x*dydt+y*dxdt
                            dstdp = x*dydp+y*dxdp
                            dsvdt = y*dzdt+z*dydt
                            dsvdp = y*dzdp+z*dydp
                            v(i, j, k) = -(cosp*dydp+sinp*dxdp) + dsvdt
                            w(i, j, k) = sinp*dzdp + cost*cosp + dstdt
                        case (2)
                            st(i, j, k) = x*z
                            sv(i, j, k) = x*y
                            dstdp = x*dzdp+z*dxdp
                            dstdt = x*dzdt+z*dxdt
                            dsvdp = x*dydp+y*dxdp
                            dsvdt = x*dydt+y*dxdt
                            !
                            !          v = -1/sin(theta)*d(st)/dphi + d(sv)/dtheta
                            !
                            !          w =  1/sin(theta)*d(sv)/dphi + d(st)/dtheta
                            !
                            v(i, j, k) = z*sinp + dsvdt
                            w(i, j, k) = cosp*dydp+ sinp*dxdp + dstdt
                    end select
                end do
            end do
        end do

        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(v(1, 1, kk), "   v", nlat, nlon)
        !     call aout(w(1, 1, kk), "   w", nlat, nlon)
        !     end do

        select case (icase)
            case (1)
        		
                call name("**ec")
        		
                call vhaeci(nlat, nlon, wsave, ierror)
                call name("vhai")
                call iout(ierror, "ierr")
        		
                call vhaec(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vha ")
                call iout(ierror, "ierr")
        		
                !     call aout(br, "  br", nlat, nlat)
                !     call aout(bi, "  bi", nlat, nlat)
                !     call aout(cr, "  cr", nlat, nlat)
                !     call aout(ci, "  ci", nlat, nlat)
        		
                !
                !     now synthesize v, w from br, bi, cr, ci and compare with original
                !
                call vhseci(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhsec(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vhs ")
                call iout(ierror, "ierr")
            case (2)
        		
                call name("**es")
        		
                call vhaesi(nlat, nlon, wsave, ierror)
                call name("vhai")
                call iout(ierror, "ierr")
        		
                call vhaes(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, ierror)
        		
                call name("vha ")
                call iout(ierror, "ierr")
        		
                !     call aout(br, "  br", nlat, nlat)
                !     call aout(bi, "  bi", nlat, nlat)
                !     call aout(cr, "  cr", nlat, nlat)
                !     call aout(ci, "  ci", nlat, nlat)
        		
                !
                !     now synthesize v, w from br, bi, cr, ci and compare with original
                !
                call vhsesi(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhses(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
        		
                call name("vhs ")
                call iout(ierror, "ierr")
            case (3)
        		
                call name("**gc")
        		
                call name("vhgi")
                call iout(nlat, "nlat")
        		
                call vhagci(nlat, nlon, wsave, ierror)
                call name("vhai")
                call iout(ierror, "ierr")
        		
                call vhagc(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vha ")
                call iout(ierror, "ierr")
        		
                !     call aout(br, "  br", nlat, nlat)
                !     call aout(bi, "  bi", nlat, nlat)
                !     call aout(cr, "  cr", nlat, nlat)
                !     call aout(ci, "  ci", nlat, nlat)
        		
                !
                !     now synthesize v, w from br, bi, cr, ci and compare with original
                !
                call vhsgci(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhsgc(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vhs ")
                call iout(ierror, "ierr")
            case (4)
        		
                call name("**gs")
        		
                call vhagsi(nlat, nlon, wsave, ierror)
                call name("vhai")
                call iout(ierror, "ierr")
        		
                call vhags(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vha ")
                call iout(ierror, "ierr")
        		
                !     call aout(br, "  br", nlat, nlat)
                !     call aout(bi, "  bi", nlat, nlat)
                !     call aout(cr, "  cr", nlat, nlat)
                !     call aout(ci, "  ci", nlat, nlat)
        		
        		
                !
                !     now synthesize v, w from br, bi, cr, ci and compare with original
                !
                call vhsgsi(nlat, nlon, wsave, ierror)
                call name("vhsi")
                call iout(ierror, "ierr")
        		
                call vhsgs(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vhs ")
                call iout(ierror, "ierr")
        end select


        !     do kk=1, nt
        !     call iout(kk, "**kk")
        !     call aout(v(1, 1, kk), "   v", nlat, nlon)
        !     call aout(w(1, 1, kk), "   w", nlat, nlon)
        !     end do

        !
        !     compute "error" in v, w
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
                    select case (k)
                        case (1)
                            st(i, j, k) = x*y
                            sv(i, j, k) = y*z
                            dstdt = x*dydt+y*dxdt
                            dstdp = x*dydp+y*dxdp
                            dsvdt = y*dzdt+z*dydt
                            dsvdp = y*dzdp+z*dydp
                            ve = -(cosp*dydp+sinp*dxdp) + dsvdt
                            we = sinp*dzdp + cost*cosp + dstdt
                        case (2)
                            st(i, j, k) = x*z
                            sv(i, j, k) = x*y
                            dstdp = x*dzdp+z*dxdp
                            dstdt = x*dzdt+z*dxdt
                            dsvdp =  x*dydp+y*dxdp
                            dsvdt = x*dydt+y*dxdt
                            ve = z*sinp + dsvdt
                            we = cosp*dydp+ sinp*dxdp + dstdt
                    end select
                    err2v = err2v + (v(i, j, k) - ve)**2
                    err2w = err2w + (w(i, j, k) - we)**2
                end do
            end do
        end do
        !
        !     set and print least squares error in v, w
        !
        err2v = sqrt(err2v/(nt*nlat*nlon))
        err2w = sqrt(err2w/(nt*nlat*nlon))
        call vout(err2v, "errv")
        call vout(err2w, "errw")
    !
    !     end of icase loop
    !
    end do

end program tvha
