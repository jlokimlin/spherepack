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
!    a program testing subroutines shpgi and shpg
!
!    shpg is the n**2 projection with complement, odd/even
!    factorization and zero truncation on a Gaussian distributed
!    grid. The projection is defined in the JCP paper "Generalized
!    discrete spherical harmonic transforms" by 
!    Paul N. Swarztrauber and William F. Spotz
!    J. Comp. Phys., 159(2000) pp. 213-230.
!
!                     April 2002
!
program tshpg

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack

    ! Explicit typing only
    implicit none

    real(wp) :: discretization_error
    real(wp) :: g
    real(wp) :: ga
    real(wp) :: gb
    real(wp) :: gh
    integer(ip) :: i
    integer(ip) :: idimg
    integer(ip) :: idp
    integer(ip) :: ierror
    integer(ip) :: iprint
    integer(ip) :: isym
    integer(ip) :: iwshp
    integer(ip) :: j
    integer(ip) :: jdimg
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
    real(wp) :: wshagc
    real(wp) :: wshp
    real(wp) :: wshsgc
    parameter (idp=8)
    parameter (kdp=idp+idp-2)
    parameter (lwshp=2*(idp+1)**2+kdp+20, &
        liwshp=4*(idp+1), lwrk=1.25*(idp+1)**2+7*idp+8)
    parameter (lwrk1=idp*kdp)
    parameter(lwork = 4*idp*(idp-1), &
        lwsha=idp*(4*idp+1)+idp+idp+15)
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp) :: work(lwrk)
    dimension sx(idp, kdp), sy(idp, kdp), &
        wshp(lwshp), iwshp(liwshp), wrk1(lwrk1)
    dimension g(idp, kdp, 2), ga(idp, idp, 2), gb(idp, idp, 2), &
        gh(idp, kdp, 2), &
        wrk2(lwork), wshagc(lwsha), wshsgc(lwsha)
    
    iprint = 0
    nt = 1
    isym = 0
    mode = 0

    do nlat=6, 8
        do mtr=1, 2
            nlon = 2*(nlat-1)
            mtrunc = nlat-mtr
            mtrunc = min(mtrunc, nlat-1, nlon/2)
            idimg = idp
            jdimg = kdp
            call shagci(nlat, nlon, wshagc, ierror)
            call check_error(ierror)

            lwshs = lwsha
            call shsgci(nlat, nlon, wshsgc, ierror)
            call check_error(ierror)

            ! Initiate faster filter
            call shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, liwshp, &
                work, lwrk, ierror)
            call check_error(ierror)

            if (iprint /= 0) write (stdout, '(3(a, i5))') &
                ' mode =' , nlon, '  nlat =', nlat, '  nlon =', nlon

            ! Initialize with pseudo random field
            do i=1, nlat
                do j=1, nlon
                    sx(i, j) =  cos(real(i*j, kind=wp))
                    g(i, j, 1) = sx(i, j)
                end do
            end do

            call shagc(nlat, nlon, mode, nt, g, idimg, jdimg, ga, gb, idimg, idimg, &
                wshagc, ierror)
            call check_error(ierror)

            if(mtrunc < nlat-1) then
                do np1=mtrunc+2, nlat
                    do mp1=1, np1
                        ga(mp1, np1, 1) = ZERO
                        gb(mp1, np1, 1) = ZERO
                    end do
                end do
            end if

            call shsgc(nlat, nlon, mode, nt, gh, idimg, jdimg, ga, gb, idimg, idimg, wshsgc, ierror)
            call check_error(ierror)

            call shpg(nlat, nlon, isym, mtrunc, sx, sy, idp, &
                wshp, lwshp, iwshp, liwshp, wrk1, lwrk1, ierror)
            call check_error(ierror)

            if (iprint > 0) write (stdout, '(/a/)') ' approx and exact solution'

            do j=1, nlon
                if(iprint > 0) write (stdout, 437) j, (sy(i, j), gh(i, j, 1), i=1, nlat)
437             format(' j=', i5, 1p4e15.6/(8x, 1p4e15.6))
            end do

            discretization_error = maxval(abs(sy(:nlat, :nlon)-gh(:nlat, :nlon, 1)))

            write (stdout, '(/2(a, i5)/)') &
                'case nlat =', nlat, ' and mtrunc =', mtrunc

            write (stdout, '(3(a, 1pe15.6/))') &
                ' error =', discretization_error, &
                ' tusl  =', tusl, &
                ' toe   =', toe
        end do
    end do
end program tshpg

