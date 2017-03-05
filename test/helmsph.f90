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
! ... file helmsph.f
!
!     this file contains a program for solving the Helmholtz
!     equation with constant 1.0 on a ten degree grid on the full sphere
!
! ... required spherepack files
!
!     islapec.f, shaec.f, shsec.f, sphcom.f, hrfft.f
!
! ... description
!
!     let theta be latitude and phi be east longitude in radians.
!     and let
!
!
!       x = cos(theta)*sin(phi)
!       y = cos(theta)*cos(phi)
!       z = sint(theta)
!
!     be the cartesian coordinates corresponding to theta and phi.
!     on the unit sphere.  The exact solution
!
!        ue(theta, phi) = (1.+x*y)*exp(z)
!
!     is used to set the right hand side and compute error.
!
!
! **********************************************************************
!
! OUTPUT FROM EXECUTING THE PROGRAM BELOW
! WITH 32 AND 64 BIT FLOATING POINT ARITHMETIC
!
! Helmholtz approximation on a ten degree grid
! nlat = 19   nlon = 36
! xlmbda =  1.00   pertrb =  0.000E+00
! maximum error =  0.715E-06 *** (32 BIT)
! maximum error =  0.114E-12 *** (64 BIT)
!
! ***********************************************
! ***********************************************
program helmsph
use spherepack
    !
    !     set grid size with parameter statements
    !
    implicit none
    integer nnlat, nnlon, nn15, llsave, llwork, lldwork
    parameter (nnlat=19, nnlon=36)
    !
    !     set saved and unsaved work space lengths in terms of nnlat, nnlon
    !     (see documentation for shaec, shsec, islapec)
    !
    parameter (nn15=nnlon+15)
    parameter (llsave=nnlat*(nnlat+1)+3*((nnlat-2)*(nnlat-1)+nn15))
    parameter (llwork=nnlat*(2*nnlon+3*(nnlat+1)+2*nnlat+1))
    !
    !     set real work space length for initializations
    !
    parameter (lldwork = nnlat+1)
    !
    !     dimension arrays
    !
    real u(nnlat, nnlon), r(nnlat, nnlon)
    real sint(nnlat), cost(nnlat), sinp(nnlon), cosp(nnlon)
    real work(llwork), wshaec(llsave), wshsec(llsave)
    real dwork(lldwork)
    real a(nnlat, nnlat), b(nnlat, nnlat)
    integer nlat, nlon, i, j, lshaec, lshsec, lwork, ierror, isym, nt
    integer ldwork
    real x, y, z, dlat, dlon, theta, phi, xlmbda(1), pertrb(1), ez, ue, errm

    !
    !     set helmholtz constant
    !
    xlmbda = 1.0
    !
    !     set work space length arguments
    !
    lwork = llwork
    ldwork = lldwork
    lshaec = llsave
    lshsec = llsave
    !
    !     set grid size arguments
    !
    nlat = nnlat
    nlon = nnlon
    !
    !     set sine and cosine vectors
    !
    dlat = pi/(nlat-1)
    dlon = (pi+pi)/nlon
    do i=1, nlat
        theta = -0.5*pi+(i-1)*dlat
        sint(i) = sin(theta)
        cost(i) = cos(theta)
    end do

    do j=1, nlon
        phi = (j-1)*dlon
        sinp(j) = sin(phi)
        cosp(j) = cos(phi)
    end do
    !
    !     set right hand side as helmholtz operator
    !     applied to ue = (1.+x*y)*exp(z)
    !
    do j=1, nlon
        do i=1, nlat
            x = cost(i)*cosp(j)
            y = cost(i)*sinp(j)
            z = sint(i)
            r(i, j) = -(x*y*(z*z+6.*(z+1.))+z*(z+2.))*exp(z)
        end do
    end do
    !
    !     initialize saved work space arrays for scalar harmonic
    !     analysis and Helmholtz inversion of r
    !
    call shaeci(nlat, nlon, wshaec, ierror)
    if (ierror > 0) then
        write (6, 200) ierror
200     format(' shaeci, ierror = ', i2)
        call exit(0)
    end if
    call shseci(nlat, nlon, wshsec, ierror)
    if (ierror > 0) then
        write (6, 201) ierror
201     format(' shseci, ierror = ', i2)
        call exit(0)
    end if
    !
    !     set no symmetry and one array
    !
    isym = 0
    nt = 1
    !
    !     compute coefficients of r for input to islapec
    !
    call shaec(nlat, nlon, isym, nt, r, nlat, nlon, a, b, nlat, nlat, &
        wshaec, ierror)
    if (ierror > 0) then
        write(*, 202) ierror
202     format(' shaec , ierror = ', i2)
        call exit(0)
    end if
    !
    !     solve Helmholtz equation on the sphere in u
    !
    write (6, 100) nlat, nlon
100 format(' helmholtz approximation on a ten degree grid' &
        /' nlat = ', i3, 2x, ' nlon = ', i3)
    call islapec(nlat, nlon, isym, nt, xlmbda, u, nlat, nlon, a, b, nlat, nlat, &
        wshsec, lshsec, work, lwork, pertrb, ierror)
    if (ierror /= 0) then
        write (6, 103) ierror
103     format(' islapec, ierror = ', i2)
        if (ierror > 0) call exit(0)
    end if
    !
    !     compute and print maximum error in u
    !
    errm = 0.0
    do j=1, nlon
        do i=1, nlat
            x = cost(i)*cosp(j)
            y = cost(i)*sinp(j)
            z = sint(i)
            ez = exp(z)
            ue = (1.+x*y)*ez
                errm = amax1(errm, abs(u(i, j)-ue))
        end do
    end do
    write(*, 204) xlmbda, pertrb, errm
204 format(' xlmbda = ', f5.2, 2x, ' pertrb = ' , e10.3, &
        /' maximum error = ', e10.3)
end program helmsph
