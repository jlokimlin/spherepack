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
! ... file testsshifte.f contains a test program illustrating the
!     use of subroutine sshifte (see documentation for sshifte)
!
! ... required spherepack files
!
!     hrfft.f, sshifte.f
!
!     The "offset" scalar field goff is set equal to exp(x+y+z) where
!     x, y, z are the cartesian coordinates restricted to the surface
!     of the unit sphere.  This is transferred to the "regular" grid
!     in greg.  greg is then compared with exp(x+y+z) on the regular grid.
!     Finally greg is transferred back to goff which is again compared
!     with exp(x+y+z).  The least squares error after each transformation
!     with sshifte is computed and printed.  Output from running the
!     program below on a 2.5 degree equally spaced regular and offset
!     grid is listed below.
!
! *** OUTPUT (from execution on 32 bit machine)
!
!           sshifte arguments
!           ioff =  0 nlon = 144 nlat =  72
!           lsave =   608 lwork = 21024
!           ier =  0
!           least squares error =  0.500E-06
!           sshifte arguments
!           ioff =  1 nlon = 144 nlat =  72
!           lsave =   608 lwork = 21024
!           ier =  0
!           least squares error =  0.666E-06
!
! *** END OF OUTPUT
!
program testsshifte
use spherepack
    implicit none
    integer nnlon, nnlat, nnlatp1, nnlat2, llsave, llwork
    parameter(nnlon=144, nnlat=72)
    parameter (nnlatp1 = nnlat+1, nnlat2 = nnlat+nnlat)
    parameter (llsave=2*(2*nnlat+nnlon)+32)
    !     for nnlon even
    parameter (llwork = 2*(nnlat+1)*nnlon)
    !     for nnlon odd
    !     parameter (llwork = nnlon*(5*nnlat+1))
    integer ioff, nlon, nlat, nlat2, j, i, lsave, lwork, ier
    real dlat, dlon, dlat2, dlon2, lat, long, x, y, z, gexact, err2
    real goff(nnlon, nnlat), greg(nnlon, nnlatp1)
    real wsave(llsave), work(llwork)
    !
    !     set resolution, work space lengths, and grid increments
    !
    nlat = nnlat
    nlon = nnlon
    lsave = llsave
    nlat2 = nnlat2
    lwork = llwork
    dlat = pi/nlat
    dlon = (pi+pi)/nlon
    dlat2 = 0.5*dlat
    dlon2 = 0.5*dlon
    !
    !     set offset grid values in goff
    !
    do j=1, nlon
        long = dlon2+(j-1)*dlon
        do i=1, nlat
            lat = -0.5*pi + dlat2+(i-1)*dlat
            x = cos(lat)*cos(long)
            y = cos(lat)*sin(long)
            z = sin(lat)
            goff(j, i) = exp(x+y+z)
        end do
    end do
    !
    !    initialize wsav for offset to regular shift
    !
    ioff = 0
    call sshifti(ioff, nlon, nlat, lsave, wsave, ier)
    !
    !     write input arguments to sshifte
    !
    write(*, 100) ioff, nlon, nlat, lsave, lwork
100 format(' sshifte arguments', &
        /' ioff = ', i2, ' nlon = ', i3, ' nlat = ', i3, &
        /' lsave = ', i5, ' lwork = ', i5)
    !
    !     shift offset to regular grid
    !
    call sshifte(ioff, nlon, nlat, goff, greg, wsave, lsave, work, lwork, ier)

    write(*, 200) ier
200 format(' ier = ', i2)

    if (ier==0) then
          !
          !     compute error in greg
          !
        err2 = 0.0
        do j=1, nlon
            long = (j-1)*dlon
            do i=1, nlat+1
                lat = -0.5*pi+(i-1)*dlat
                x = cos(lat)*cos(long)
                y = cos(lat)*sin(long)
                z = sin(lat)
                gexact = exp(x+y+z)
                err2 = err2 + (greg(j, i)-gexact)**2
            end do
        end do
        err2 = sqrt(err2/(nlon*(nlat+1)))
        write(*, 300) err2
300     format(' least squares error = ', e10.3)
    end if
    !    initialize wsav for regular to offset shift
    !
    ioff = 1
    call sshifti(ioff, nlon, nlat, lsave, wsave, ier)
    !
    !     now transfer regular grid values in greg back to offset grid in goff
    !
    do j=1, nlon
        do i=1, nlat
            goff(j, i) = 0.0
        end do
    end do
    write(*, 100) ioff, nlon, nlat, lsave, lwork
    call sshifte(ioff, nlon, nlat, goff, greg, wsave, lsave, work, lwork, ier)
    write(*, 200) ier
    if (ier == 0) then
        !
        !     compute error in goff by comparing with exp(x+y+z) on offset grid
        !
        err2 = 0.0
        do j=1, nlon
            long = dlon2+(j-1)*dlon
            do i=1, nlat
                lat = -0.5*pi+dlat2+(i-1)*dlat
                x = cos(lat)*cos(long)
                y = cos(lat)*sin(long)
                z = sin(lat)
                gexact = exp(x+y+z)
                err2 = err2 + (goff(j, i)-gexact)**2
            end do
        end do
        err2 = sqrt(err2/(nlon*(nlat+1)))
        write(*, 300) err2
    end if
end program testsshifte
