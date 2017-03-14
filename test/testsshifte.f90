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
! ... Contains a test program illustrating the
!     use of subroutine sshifte (see documentation for sshifte)
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

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack

    ! Explicit typing only
    implicit none

    integer(ip), parameter :: nlon = 144, nlat = 72
    integer(ip), parameter :: nlatp1 = nlat + 1
    integer(ip) :: ioff, i, j, error_flag
    real(wp) :: dlat, dlon, half_dlat, half_dlon, lat, long
    real(wp) :: x, y, z, gexact, err2
    real(wp) :: goff(nlon, nlat), greg(nlon, nlatp1)
    real(wp), allocatable :: wavetable(:)
    real(wp), parameter :: ZERO = 0.0_wp

    ! Set grid increments
    dlat = PI/nlat
    dlon = TWO_PI/nlon
    half_dlat = dlat/2
    half_dlon = dlon/2

    ! Set offset grid values in goff
    do j=1, nlon
        long = half_dlon + real(j - 1, kind=wp) * dlon
        do i=1, nlat
            lat = -HALF_PI + half_dlat + real(i - 1, kind=wp) * dlat
            x = cos(lat) * cos(long)
            y = cos(lat) * sin(long)
            z = sin(lat)
            goff(j, i) = exp(x + y + z)
        end do
    end do

    ! Initialize wavetable for offset to regular shift
    ioff = 0
    call initialize_sshifte(ioff, nlon, nlat, wavetable, error_flag)

    ! Write input arguments to sshifte
    write (stdout, 100) ioff, nlon, nlat
100 format(' sshifte arguments', &
        /' ioff = ', i2, ' nlon = ', i3, ' nlat = ', i3)

    ! Shift offset to regular grid
    call sshifte(ioff, nlon, nlat, goff, greg, wavetable, error_flag)

    write (stdout, 200) error_flag
200 format(' ier = ', i2)

    if (error_flag == 0) then
        !
        !     compute error in greg
        !
        err2 = ZERO
        do j=1, nlon
            long = real(j - 1, kind=wp) * dlon
            do i=1, nlat + 1
                lat = -HALF_PI + real(i - 1, kind=wp) * dlat
                x = cos(lat)*cos(long)
                y = cos(lat)*sin(long)
                z = sin(lat)
                gexact = exp(x + y + z)
                err2 = err2 + (greg(j, i)-gexact)**2
            end do
        end do
        err2 = sqrt(err2/(nlon*(nlat + 1)))
        write (stdout, 300) err2
300     format(' least squares error = ', e10.3)
    end if

    ! Initialize wsav for regular to offset shift
    ioff = 1
    call initialize_sshifte(ioff, nlon, nlat, wavetable, error_flag)

    ! Now transfer regular grid values in greg back to offset grid in goff
    goff = ZERO

    write (stdout, 100) ioff, nlon, nlat

    call sshifte(ioff, nlon, nlat, goff, greg, wavetable, error_flag)

    write (stdout, 200) error_flag

    if (error_flag == 0) then

        ! Compute error in goff by comparing with exp(x+y+z) on offset grid
        err2 = ZERO
        do j=1, nlon
            long = half_dlon + real(j - 1, kind=wp) * dlon
            do i=1, nlat
                lat = -HALF_PI + half_dlat + real(i - 1, kind=wp) * dlat
                x = cos(lat) * cos(long)
                y = cos(lat) * sin(long)
                z = sin(lat)
                gexact = exp(x + y + z)
                err2 = err2 + (goff(j, i) - gexact)**2
            end do
        end do
        err2 = sqrt(err2/(nlon*(nlat + 1)))
        write (stdout, 300) err2
    end if

    ! Release memory
    deallocate (wavetable)

end program testsshifte
