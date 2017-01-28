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
!    a program testing subroutines shpei and shpe in file shpe.f
!    also requires files  shaec.f shsec.f sphcom.f and hrfft.f
!
!    shpe is the n**2 filter with complement, odd/even
!    factorization and zero truncation on an equally spaced
!    grid. The projection is defined in the JCP paper "Generalized
!    discrete spherical harmonic transforms" by 
!    Paul N. Swarztrauber and William F. Spotz
!    J. Comp. Phys., 159(2000) pp. 213-230.
!
!                     April 2002
!
program tshpe

    use, intrinsic :: ISO_Fortran_env, only: &
        sp => REAL32, &
        stderr => ERROR_UNIT, &
        stdout => OUTPUT_UNIT

    use spherepack_library

    ! Explicit typing only
    implicit none

    real(wp) :: discretization_error
    real(wp) :: g
    real(wp) :: ga
    real(wp) :: gb
    real(wp) :: gh
    real(wp) :: gw
    integer(ip) :: i
    integer(ip) :: idimg
    integer(ip) :: idp
    integer(ip) :: ierror
    integer(ip) :: iprint
    integer(ip) :: isym
    integer(ip) :: iwshp
    integer(ip) :: j
    integer(ip) :: jdimg
    integer(ip) :: k
    integer(ip) :: kdp
    integer(ip) :: liwshp
    integer(ip) :: lwork
    integer(ip) :: lwrk
    integer(ip) :: lwrk1
    integer(ip) :: lwsha
    integer(ip) :: lwshp
    integer(ip) :: lwshs
    integer(ip) :: mode
    integer(ip) :: mp1
    integer(ip) :: mtr
    integer(ip) :: mtrunc
    integer(ip) :: nlat
    integer(ip) :: nlon
    integer(ip) :: np1
    integer(ip) :: nt
    real(wp) :: sx
    real(wp) :: sy
    
    real(wp) :: toe
    real(wp) :: tusl
    real(wp) :: wrk1
    real(wp) :: wrk2
    real(wp) :: wshaec
    real(wp) :: wshp
    real(wp) :: wshsec
    parameter (idp=32)
    parameter (kdp=idp+idp-2)
    parameter (lwshp=2*(idp+1)**2+kdp+20, &
        liwshp=4*(idp+1), lwrk=1.25*(idp+1)**2+7*idp+8)
    parameter (lwrk1=idp*kdp)
    parameter(lwork = 5*idp*(idp-1), &
        lwsha=idp*(idp+1)+3*(idp-2)*(idp-1)/2+kdp+15)
    real(wp) :: work(lwrk)
    dimension sx(idp, kdp), sy(idp, kdp), &
        wshp(lwshp), iwshp(liwshp), wrk1(lwrk1)
    dimension g(idp, kdp, 2), ga(idp, idp, 2), gb(idp, idp, 2), &
        gh(idp, kdp, 2), gw(idp, kdp, 2), &
        wrk2(lwork), wshaec(lwsha), wshsec(lwsha)
    real(wp), parameter :: ZERO = 0.0_wp

    iprint = 0
    nt = 1

    do nlat=6, 8
        do mtr=1, 3
            mtrunc = nlat-mtr
            idimg = idp
            jdimg = kdp
            nlon = 2*(nlat-1)
            mtrunc = min(mtrunc, nlat-1, nlon/2)

            call shaeci(nlat, nlon, wshaec, lwsha, work, lwrk, ierror)
            call check_error(ierror)

            lwshs = lwsha

            call shseci(nlat, nlon, wshsec, lwshs, work, lwrk, ierror)
            call check_error(ierror)

            mode = 0
            if (iprint /= 0) write(stdout, '(3(a, i5))') &
                ' mode =' , nlon, '  nlat =', nlat, '  nlon =', nlon

            do k=1, nt
                do i=1, nlat
                    do j=1, nlon
                        g(i, j, k) = cos(real(i*j, kind=wp))
                        gh(i, j, k) = g(i, j, k)
                    end do
                end do
            end do

            call shaec(nlat, nlon, mode, nt, g, idimg, jdimg, ga, gb, idimg, idimg, &
                wshaec, lwsha, wrk2, lwork, ierror)
            call check_error(ierror)

            if (mtrunc < nlat-1) then
                do np1=mtrunc+2, nlat
                    do mp1=1, np1
                        ga(mp1, np1, 1) = ZERO
                        gb(mp1, np1, 1) = ZERO
                    end do
                end do
            end if

            call shsec(nlat, nlon, mode, nt, gw, idimg, jdimg, ga, gb, idimg, idimg, &
                wshsec, lwshs, wrk2, lwork, ierror)
            call check_error(ierror)

            ! faster filter

            isym = 0
            call shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, liwshp, &
                work, lwrk, ierror)
            call check_error(ierror)

            do j=1, nlon
                do i=1, nlat
                    sx(i, j) = gh(i, j, 1)
                end do
            end do

            call shpe(nlat, nlon, isym, mtrunc, sx, sy, idp, &
                wshp, lwshp, iwshp, liwshp, wrk1, lwrk1, ierror)
            call check_error(ierror)

            if (iprint > 0) write(stdout, '(/a/)') ' approx and exact solution'

            do j=1, nlon
                if (iprint > 0) write(stdout, 437) j, (sy(i, j), gw(i, j, 1), i=1, nlat)
437             format(' j=', i5, 1p4e15.6/(8x, 1p4e15.6))
            end do

            discretization_error = maxval(abs(sy(:nlat,:nlon)-gw(:nlat,:nlon,1)))

            write(stdout, '(/2(a, i5)/)') &
                'case nlat =', nlat, ' and mtrunc =', mtrunc

            write(stdout, '(3(a,1pe15.6/))') &
                ' error =', discretization_error, &
                ' tusl  =', tusl, &
                ' toe   =', toe
        end do
    end do

end program tshpe

