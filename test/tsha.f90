!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                          Spherepack                           *
!     *                                                               *
!     *       A Package of Fortran subroutines and programs           *
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
!     a program for testing all scalar analysis and synthesis subroutines
!
program test_all_scalar_analysis_and_synthesis_routines

    use spherepack

    implicit none

    !
    !     set dimensions with parameter statements
    !
    integer(ip), parameter :: nlat= 15, nlon= 18, nt = 3
    integer(ip), parameter :: lleng= 5 * nlat * nlat * nlon, llsav = 5*nlat*nlat*nlon
    integer(ip), parameter :: lldwork = nlat * (nlat + 4)

    real(wp) :: cosp
    real(wp) :: cost
    real(wp) :: dlat
    real(wp) :: dphi
    real(wp) :: err2
    integer(ip) :: i
    integer(ip) :: icase
    integer(ip) :: error_flag
    integer(ip) :: ierror
    integer(ip) :: isym
    integer(ip) :: j
    integer(ip) :: k
    
    integer(ip) :: ldwork
    integer(ip) :: lsave
    integer(ip) :: lwork
    
    real(wp) :: phi
    real(wp) :: sinp
    real(wp) :: sint
    real(wp) :: theta
    real(wp) :: xyzk

    real(wp) :: dwork(lldwork), work(lleng), wsave(llsav)
    real(wp)                          :: s(nlat,nlon,nt)
    real(wp), dimension(nlat,nlat,nt) :: a, b
    real(wp), dimension(nlat)         :: gaussian_latitudes, gaussian_weights
    !
    !     set dimension variables
    !
    lwork = lleng
    lsave = llsav
    call iout(nlat, "nlat")
    call iout(nlon, "nlon")
    call iout(nt, "  nt")
    isym = 0
    !
    !     set equally spaced colatitude and longitude increments
    !
    dphi = TWO_PI/nlon
    dlat = pi/(nlat-1)
    !
    !     compute nlat gaussian points in thetag
    !
    ldwork = lldwork
    call compute_gaussian_latitudes_and_weights(nlat, gaussian_latitudes, gaussian_weights, error_flag)

    call name("gaqd")
    call iout(error_flag, " ier")
    call vecout(gaussian_latitudes, "thtg", nlat)
    !
    !     test all analysis and synthesis subroutines
    !
    do icase=1, 4
        !
        !     icase=1 test shaec, shsec
        !     icase=2 test shaes, shses
        !     icase=3 test shagc, shsgc
        !     icase=4 test shags, shsgs
        !
        call name("****")
        call name("****")
        call iout(icase, "icas")
        !
        !
        !     set scalar field as (x*y*z)**k) restricted to the sphere
        !
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    xyzk = (sint*(sint*cost*sinp*cosp))**k
                    !           s(i, j, k) = exp(xyzk)
                    s(i, j, k) = xyzk
                end do
            end do
        !     call iout(k, "   k")
        !     call aout(s(1, 1, k), "   s", nlat, nlon)
        end do

        wsave(1: lsave) = 0.0

        select case (icase)
            case (1)
        		
                call name("**ec")
                call shaeci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shaec(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("sha ")
                call iout(ierror, "ierr")
        		
                call shseci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		
                call name("shsi")
                call iout(ierror, "ierr")
        		
                call shsec(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("shs ")
                call iout(ierror, "ierr")
            case (2)
        		
                call name("**es")
                call shaesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shaes(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("sha ")
                call iout(ierror, "ierr")
        		
                call shsesi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
        		
                call name("shsi")
                call iout(ierror, "ierr")
        		
                call shses(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("shs ")
                call iout(ierror, "ierr")
            case (3)
        		
                call name("**gc")
        		
                call shagci(nlat, nlon, wsave, ierror)
        		
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shagc(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("sha ")
                call iout(ierror, "ierr")
        		
                call shsgci(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
        		
                call name("shsi")
                call iout(ierror, "ierr")
        		
                call shsgc(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("shs ")
                call iout(ierror, "ierr")
            case (4)
        		
                call name("**gs")
        		
                call shagsi(nlat, nlon, wsave, ierror)
        		
                call name("shai")
                call iout(ierror, "ierr")
        		
                call shags(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("sha ")
                call iout(ierror, "ierr")
        		
                call shsgsi(nlat, nlon, wsave, lsave, work, lwork, dwork, ldwork, ierror)
                call name("shsi")
                call iout(ierror, "ierr")
        		
                call shsgs(nlat, nlon, isym, nt, s, nlat, nlon, a, b, nlat, nlat, wsave, &
                    lsave, work, lwork, ierror)
        		
                call name("shs ")
                call iout(ierror, "ierr")
        end select
        !
        !     compute "error" in s
        !
        err2 = 0.0
        do k=1, nt
            do j=1, nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1, nlat
                    theta = (i-1)*dlat
                    if (icase > 2) theta = gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    xyzk = (sint*(sint*cost*sinp*cosp))**k
                    !           err2 = err2+ (exp(xyzk)-s(i, j, k))**2
                    err2 = err2 + (xyzk-s(i, j, k))**2
                end do
            end do
        !     call iout(k, "   k")
        !     call aout(s(1, 1, k), "   s", nlat, nlon)
        end do
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2, "err2")
    end do

end program test_all_scalar_analysis_and_synthesis_routines
