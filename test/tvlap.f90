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
!     1/976
!
!     a program for testing all vector laplacian and its inverse
!
!     (1) first set the vector function (rotation about equator) using
!         v = cos(phi) and w = -cos(theta)*sin(phi)
!
!     (2) set vector laplacian ananlytically
!         vlap = -2.*cos(phi)=-2.*v, wlap = -2.*w
!         (i.e., L(v, w) = -2.*(v, w) so (v, w) is an eigenfunction for the
!         vector Laplacian with eigenvalue -2.
!
!     (3) compute the coefficients br, bi, cr, ci of (v, w) using vector analysis
!
!     (3) compute the vector laplacian of (v, w) using vlapec, vlapes, vlapgc, vlapgs
!
!     (4) compare (3) with (2)
!
!     (5) invert (4) and compare with (v, w)
!
program tvlap
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
    integer(ip) :: mdbc
    integer(ip) :: mmdbc
    integer(ip) :: nlat
    integer(ip) :: nlon
    integer(ip) :: nmax
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
    real(wp) :: ve
    real(wp) :: vlap
    real(wp) :: w
    real(wp) :: we
    real(wp) :: wlap
    real(wp) :: work
    real(wp) :: wsave
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat=29 , nnlon= 16, nnt = 1)
    parameter (mmdbc = (nnlon+2)/2)
    parameter (lleng= 5*nnlat*nnlat*nnlon, llsav=15*nnlat*nnlat*nnlon)
    parameter (lldwork = 4*nnlat*nnlat)
    
    dimension work(lleng), wsave(llsav)
    dimension br(mmdbc, nnlat, nnt), bi(mmdbc, nnlat, nnt)
    dimension cr(mmdbc, nnlat, nnt), ci(mmdbc, nnlat, nnt)
    dimension thetag(nnlat), dtheta(nnlat), dwts(nnlat)
    dimension v(nnlat, nnlon, nnt), w(nnlat, nnlon, nnt)
    dimension vlap(nnlat, nnlon, nnt), wlap(nnlat, nnlon, nnt)

    real dtheta, dwts
    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    nmax = max(nlat, nlon)
    mdbc = mmdbc

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
                    if (k==1) then
                        v(i, j, k) = cosp
                        w(i, j, k) = -cost*sinp
                        vlap(i, j, k) = -2.0*v(i, j, k)
                        wlap(i, j, k) = -2.0*w(i, j, k)
                    end if
                end do
            end do
        end do

        if (nmax<10) then
            do kk=1, nt
                call iout(kk, "**kk")
            !     call aout(v(1, 1, kk), "   v", nlat, nlon)
            !     call aout(w(1, 1, kk), "   w", nlat, nlon)
            !     call aout(vlap(1, 1, kk), "vlap", nlat, nlon)
            !     call aout(wlap(1, 1, kk), "wlap", nlat, nlon)
            end do
        end if

        if (icase==1) then

            call name("**ec")
            !
            !     analyze vector field
            !
            call vhaeci(nlat, nlon, wsave, ierror)
            call name("vhai")
            call iout(ierror, "ierr")

            call vhaec(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                nlat, wsave, ierror)
            call name("vha ")
            call iout(ierror, "ierr")

            !     if (nmax.lt.10) then
            !     do kk=1, nt
            !     call iout(kk, "**kk")
            !     call aout(br(1, 1, kk), "  br", nlat, nlat)
            !     call aout(bi(1, 1, kk), "  bi", nlat, nlat)
            !     call aout(cr(1, 1, kk), "  cr", nlat, nlat)
            !     call aout(ci(1, 1, kk), "  ci", nlat, nlat)
            !     end do
            !     end if
            !
            !     compute vector laplacian
            !

            call vhseci(nlat, nlon, wsave, ierror)
            call name("vhsi")
            call iout(ierror, "ierr")

            call vlapec(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("vlap")
            call iout(ierror, "ierr")


        else if (icase==2) then

            call name("**es")
            !
            !     analyze vector field
            !
            call vhaesi(nlat, nlon, wsave, ierror)
            call name("vhaesi")
            call iout(ierror, "ierr")

            call vhaes(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                nlat, wsave, ierror)
            call name("vhaes")
            call iout(ierror, "ierr")

            call vhsesi(nlat, nlon, wsave, ierror)
            call name("vhsesi")
            call iout(ierror, "ierr")

            call vlapes(nlat, nlon, isym, nt, vlap, wlap, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("vlapes")
            call iout(ierror, "ierr")

        else if (icase ==3) then

            call name("**gc")
            !
            !     analyze vector field
            !
            call vhagci(nlat, nlon, wsave, ierror)
            call name("vhagci")
            call iout(ierror, "ierr")

            call vhagc(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vhagc")
            call iout(ierror, "ierr")

            !     if (nmax.lt.10) then
            !     do kk=1, nt
            !     call iout(kk, "**kk")
            !     call aout(br(1, 1, kk), "  br", nlat, nlat)
            !     call aout(bi(1, 1, kk), "  bi", nlat, nlat)
            !     call aout(cr(1, 1, kk), "  cr", nlat, nlat)
            !     call aout(ci(1, 1, kk), "  ci", nlat, nlat)
            !     end do
            !     end if
            !
            !     compute vector laplacian
            !

            call vhsgci(nlat, nlon, wsave, ierror)
            call name("vhsgci")
            call iout(ierror, "ierr")

            call vlapgc(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("vlapgc")
            call iout(ierror, "ierr")

        else if (icase == 4) then

            call name("**gs")
            !
            !     analyze vector field
            !
            call vhagsi(nlat, nlon, wsave, ierror)
            call name("vhagsi")
            call iout(ierror, "ierr")

            call vhags(nlat, nlon, ityp, nt, v, w, nlat, nlon, br, bi, cr, ci, mdbc, &
                nlat, wsave, lsave, work, lwork, ierror)
            call name("vhags")
            call iout(ierror, "ierr")

            !     if (nmax.lt.10) then
            !     do kk=1, nt
            !     call iout(kk, "**kk")
            !     call aout(br(1, 1, kk), "  br", nlat, nlat)
            !     call aout(bi(1, 1, kk), "  bi", nlat, nlat)
            !     call aout(cr(1, 1, kk), "  cr", nlat, nlat)
            !     call aout(ci(1, 1, kk), "  ci", nlat, nlat)
            !     end do
            !     end if
            !
            !     compute vector laplacian
            !

            call vhsgsi(nlat, nlon, wsave, ierror)
            call name("vhsgsi")
            call iout(ierror, "ierr")

            call vlapgs(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("vlapgs")
            call iout(ierror, "ierr")

        end if

        if (nmax<10) then
            do kk=1, nt
                call iout(kk, "**kk")
            !     call aout(vlap(1, 1, kk), "vlap", nlat, nlon)
            !     call aout(wlap(1, 1, kk), "wlap", nlat, nlon)
            end do
        end if
        !
        !     compute "error" in vlap, wlap
        !
        err2v = 0.0
        err2w =0.0
        do k=1, nt
            do j=1, nlon
                do i=1, nlat
                    if (k==1) then
                        err2v = err2v+(vlap(i, j, k)+2.*v(i, j, k))**2
                        err2w = err2w+(wlap(i, j, k)+2.*w(i, j, k))**2
                    end if
                end do
            end do
        end do
        !
        !     set and print least squares error in vlap, wlap
        !
        err2v = sqrt(err2v/(nt*nlat*nlon))
        err2w = sqrt(err2w/(nt*nlon*nlat))
        call vout(err2v, "errv")
        call vout(err2w, "errw")
        !
        !     now recompute (v, w) inverting (vlap, wlap) ivlap codes
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

            call name("analyze vector field (vlap, wlap)")
            !
            !     analyze vector field (vlap, wlap)
            !
            call vhaeci(nlat, nlon, wsave, ierror)
            call name("vhaeci")
            call iout(ierror, "ierr")

            call vhaec(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, &
                br, bi, cr, ci, mdbc, nlat, wsave, ierror)
            call name("vhaec ")
            call iout(ierror, "ierr")

            call vhseci(nlat, nlon, wsave, ierror)
            call name("vhseci")
            call iout(ierror, "ierr")

            call ivlapec(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("ivlapec")
            call iout(ierror, "ierr")

        else if (icase==2) then

            call name("analyze vector field (vlap, wlap)")
            !
            !     analyze vector field (vlap, wlap)
            !
            call vhaesi(nlat, nlon, wsave, ierror)
            call name("vhaesi")
            call iout(ierror, "ierr")

            call vhaes(nlat, nlon, isym, nt, vlap, wlap, nlat, nlon, &
                br, bi, cr, ci, mdbc, nlat, wsave, ierror)
            call name("vhaes")
            call iout(ierror, "ierr")

            call vhsesi(nlat, nlon, wsave, ierror)
            call name("vhsesi")
            call iout(ierror, "ierr")

            call ivlapes(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("ivlapes")
            call iout(ierror, "ierr")

        else if (icase == 3) then

            call name("analyze vector field (vlap, wlap)")

            !
            !     analyze vector field (vlap, wlap)
            !
            call vhagci(nlat, nlon, wsave, ierror)
            call name("vhagci")
            call iout(ierror, "ierr")

            call vhagc(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, &
                br, bi, cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("vhagc")
            call iout(ierror, "ierr")

            call vhsgci(nlat, nlon, wsave, ierror)
            call name("vhsgci")
            call iout(ierror, "ierr")

            call ivlapgc(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("ivlapgc")
            call iout(ierror, "ierr")

        else if (icase == 4) then

            call name("analyze vector field (vlap, wlap)")

            !
            !     analyze vector field (vlap, wlap)
            !
            call vhagsi(nlat, nlon, wsave, ierror)
            call name("vhagsi")
            call iout(ierror, "ierr")

            call vhags(nlat, nlon, ityp, nt, vlap, wlap, nlat, nlon, &
                br, bi, cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("vhags")
            call iout(ierror, "ierr")

            call vhsgsi(nlat, nlon, wsave, ierror)
            call name("vhsgsi")
            call iout(ierror, "ierr")

            call ivlapgs(nlat, nlon, isym, nt, v, w, nlat, nlon, br, bi, &
                cr, ci, mdbc, nlat, wsave, lsave, work, lwork, ierror)
            call name("ivlapgs")
            call iout(ierror, "ierr")

        end if


        if (nmax<10) then
            do kk=1, nt
                call iout(kk, "**kk")
            !     call aout(v(1, 1, kk), "   v", nlat, nlon)
            !     call aout(w(1, 1, kk), "   w", nlat, nlon)
            end do
        end if

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
                    if (k==1) then
                        ve = cosp
                        we = -cost*sinp
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

end program tvlap
