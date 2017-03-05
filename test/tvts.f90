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
!     4/97
!
!     a program for testing all theta derivative subroutines
!     vtses, vtsec, vtsgs, vtsgc
!
!
!     (1) first set a valid vector field (v, w) in terms of x, y, z
!         cartesian coordinates
!
!     (2) analytically compute (vt, wt) from (1)
!
!     (3) compute (vt, wt) using vtses, vtsec, vtsgs, vtsgc and compare with (2)
!
program tvts
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
    real(wp) :: dxdt
    real(wp) :: dydt
    real(wp) :: dzdt
    real(wp) :: emz
    real(wp) :: err2v
    real(wp) :: err2w
    real(wp) :: ex
    real(wp) :: ey
    real(wp) :: ez
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
    real(wp) :: theta
    real(wp) :: thetag
    real(wp) :: v
    real(wp) :: vt
    real(wp) :: vtsav
    real(wp) :: w
    real(wp) :: work
    real(wp) :: wsave
    real(wp) :: wt
    real(wp) :: wtsav
    real(wp) :: x
    real(wp) :: y
    real(wp) :: z
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 25, nnlon= 19, nnt = 3)
    parameter (lleng= 5*nnlat*nnlat*nnlon, llsav= 5*nnlat*nnlat*nnlon)
    dimension work(lleng), wsave(llsav)
    parameter (lldwork = 4*nnlat*nnlat)
    real dwork(lldwork)
    dimension br(nnlat, nnlat, nnt), bi(nnlat, nnlat, nnt)
    dimension cr(nnlat, nnlat, nnt), ci(nnlat, nnlat, nnt)
    dimension thetag(nnlat), dtheta(nnlat), dwts(nnlat)
    dimension v(nnlat, nnlon, nnt), w(nnlat, nnlon, nnt)
    dimension vt(nnlat, nnlon, nnt), wt(nnlat, nnlon, nnt)
    dimension vtsav(nnlat, nnlon, nnt), wtsav(nnlat, nnlon, nnt)
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
    dphi = TWO_PI/nlon
    dlat = PI/(nlat-1)
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
    !     test all theta derivative subroutines
    !
    do icase=1, 4
        !
        !     icase=1 test vtsec
        !     icase=2 test vtses
        !     icase=3 test vtsgc
        !     icase=4 test vtsgs
        !
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !
        !     set vector field v, w and compute theta derivatives in (vtsav, wtsav)
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
                       !
                       !    set x, y, z and their theta derivatives at colatitude theta and longitude p
                       !
                    x = sint*cosp
                    dxdt = cost*cosp
                    y = sint*sinp
                    dydt = cost*sinp
                    z = cost
                    dzdt = -sint
                    !
                    !     set (v, w) field corresponding to stream function
                    !     S = exp(y)+exp(-z) and velocity potential function
                    !     P = exp(x)+exp(z)
                    !
                    ex = exp(x)
                    ey = exp(y)
                    ez = exp(z)
                    emz = exp(-z)
                    w(i, j, k) =-ex*sinp+emz*sint+ey*cost*sinp
                    v(i, j, k) =-ey*cosp-ez*sint+ex*cost*cosp
                    !
                    !     set theta derivatives differentiating w, v above
                    !
                    wtsav(i, j, k) = -ex*dxdt*sinp+emz*(-dzdt*sint+cost) &
                        +ey*sinp*(dydt*cost-sint)
                    vtsav(i, j, k) = -ey*dydt*cosp-ez*(dzdt*sint+cost) &
                        +ex*cosp*(dxdt*cost-sint)
                end do
            end do
        end do

        !     call a3out(wtsav, "wtsv", nlat, nlon, nt)
        !     call a3out(vtsav, "vtsv", nlat, nlon, nt)



        select case (icase)
            case (1)
        		
                call name("**ec")
        		
                call vhaeci(nlat, nlon, wsave, ierror)
                call name("vhai")
                call iout(ierror, "ierr")
        		
                call vhaec(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, ierror)
                call name("vha ")
                call iout(ierror, "ierr")
        		
                !
                !     now compute theta derivatives of v, w
                !
                call vtseci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		
                call name("vtsi")
                call iout(ierror, "ierr")
        		
                call vtsec(nlat, nlon, ityp, nt, vt, wt, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vts ")
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
        		
                call vtsesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
                call name("vtsi")
                call iout(ierror, "ierr")
        		
                call vtses(nlat, nlon, ityp, nt, vt, wt, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vts ")
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
        		
                !
                !     now synthesize v, w from br, bi, cr, ci and compare with original
                !
                call vtsgci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
                call name("vtsi")
                call iout(ierror, "ierr")
        		
                call vtsgc(nlat, nlon, ityp, nt, vt, wt, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vts ")
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
        		
                call vtsgsi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
                call name("vtsi")
                call iout(ierror, "ierr")
        		
                call vtsgs(nlat, nlon, ityp, nt, vt, wt, nlat, nlon, br, bi, cr, ci, nlat, &
                    nlat, wsave, lsave, work, lwork, ierror)
                call name("vts ")
                call iout(ierror, "ierr")
        end select

        !     call a3out(wt, "wt  ", nlat, nlon, nt)
        !     call a3out(vt, "vt  ", nlat, nlon, nt)

        !
        !     compute "error" in vt, wt
        !
        err2v = 0.0
        err2w = 0.0
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    err2v = err2v + (vt(i, j, k) - vtsav(i, j, k))**2
                    err2w = err2w + (wt(i, j, k) - wtsav(i, j, k))**2
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

end program tvts
