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
!
! ... A test program illustrating the
!     use of subroutine vshifte (see documentation for vshifte)
!
!
!     Let the analytic vector field (u, v) in geophysical coordinates be
!     given by
!
!          u  = -exp(x)*sin(p) + exp(-z)*cos(p) + exp(y)*sin(t)*sin(p)
!
!          v  = exp(y)*cos(p) + exp(z)*cos(t) - exp(x)*sin(t)*cos(p)
!
!     where t is the latitude coordinate, p is the longitude coordinate, 
!     and x=cos(t)*cos(p), y=cos(t)*sin(p), z=sin(t) are the cartesian coordinates
!     restricted to the sphere.

!     The "offset" vector field (uoff, voff) is set equal to (u, v).
!     This is transferred to the "regular" grid in (ureg, vreg).  (ureg, vreg)
!     is then compared with (u, v) on the regular grid.  Finally (ureg, vreg)
!     is transferred back to (uoff, voff) which is again compared with (u, v).
!     The least squares error after each transformation with vshifte
!     is computed and printed.  Results from running the program on
!     a 2.5 degree equally spaced regular and offset grid is given.
!     Output from runs on separate platforms with 32 bit and 64 bit
!     floating point arithmetic is listed.
!
! *********************************************
!     OUTPUT
! *********************************************
!
!       vshifte arguments
!       ioff =  0 nlon = 144 nlat =  72
!       lsave =   608 lwork = 21024
!       ier =  0
!       least squares error
! ***   32 BIT ARITHMETIC
!       err2u =  0.377E-06 err2v =  0.328E-06
! ***   64 BIT ARITHMETIC
!       err2u =  0.777E-13 err2v =  0.659E-13

!       vshifte arguments
!       ioff =  1 nlon = 144 nlat =  72
!       lsave =   608 lwork = 21024
!       ier =  0
!       least squares error
! ***   32 BIT ARITHMETIC
!       err2u =  0.557E-06 err2v =  0.434E-06
! ***   64 BIT AIRTHMETIC
!       err2u =  0.148E-12 err2v =  0.118E-12
!
! *********************************************
!     END OF OUTPUT (CODE FOLLOWS)
! *********************************************
!
program testvshifte

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack

    ! Explicit typing only
    implicit none

    ! Dictionary
    integer(ip), parameter :: nlon = 144, nlat = 72
    integer(ip), parameter :: nlatp1 = nlat+1
    integer(ip) :: ioff, i, j, error_flag
    real(wp) :: dlat, dlon, half_dlat, half_dlon, lat, long, x, y, z, ex, ey, ez, emz
    real(wp) :: err2u, err2v, ue, ve, sint, sinp, cost, cosp
    real(wp) :: uoff(nlon, nlat), voff(nlon, nlat)
    real(wp) :: ureg(nlon, nlatp1), vreg(nlon, nlatp1)
    real(wp), allocatable :: wavetable(:)
    real(wp), parameter   :: ZERO = 0.0_wp

    ! Set grid increments
    dlat = PI/nlat
    dlon = TWO_PI/nlon
    half_dlat = dlat/2
    half_dlon = dlon/2

    ! Set (uoff, voff) = (u, v) on offset grid
    do j=1, nlon
        long = half_dlon + real(j - 1, kind=wp) * dlon
        sinp = sin(long)
        cosp = cos(long)
        do i=1, nlat
            lat = -HALF_PI + half_dlat + real(i - 1, kind=wp) * dlat
            sint = sin(lat)
            cost = cos(lat)
            x = cost*cosp
            y = cost*sinp
            z = sint
            ex = exp(x)
            ey = exp(y)
            ez = exp(z)
            emz = exp(-z)
            uoff(j, i) =-ex*sinp+emz*cost+ey*sint*sinp
            voff(j, i) = ey*cosp+ez*cost-ex*sint*cosp
        end do
    end do

    ! Initialize wsav for offset to regular shift
    ioff = 0
    call initialize_vshifte(ioff, nlon, nlat, wavetable, error_flag)

    ! Write input arguments to vshifte
    write (stdout, 100) ioff, nlon, nlat
100 format(' vshifte arguments', &
        /' ioff = ', i2, ' nlon = ', i3, ' nlat = ', i3)

    ! Shift offset to regular grid
    call vshifte(ioff, nlon, nlat, uoff, voff, ureg, vreg, wavetable, error_flag)

    write (stdout, 200) error_flag
200 format(' ier = ', i2)

    if (error_flag == 0) then
        ! Compute error in ureg, vreg
        err2u = ZERO
        err2v = ZERO
        do j=1, nlon
            long = real(j - 1, kind=wp) * dlon
            sinp = sin(long)
            cosp = cos(long)
            do i=1, nlat+1
                lat = -HALF_PI + real(i - 1, kind=wp) * dlat
                sint = sin(lat)
                cost = cos(lat)
                x = cost*cosp
                y = cost*sinp
                z = sint
                ex = exp(x)
                ey = exp(y)
                ez = exp(z)
                emz = exp(-z)
                ue = -ex*sinp+emz*cost+ey*sint*sinp
                ve = ey*cosp+ez*cost-ex*sint*cosp
                err2u = err2u + (ureg(j, i)-ue)**2
                err2v = err2v + (vreg(j, i)-ve)**2
            end do
        end do
        err2u = sqrt(err2u/(nlon*(nlat + 1)))
        err2v = sqrt(err2v/(nlon*(nlat + 1)))
        write (stdout, 300) err2u, err2v
300     format(' least squares error ', &
            /' err2u = ', e10.3, ' err2v = ', e10.3)
    end if

    ! Initialize wsav for regular to offset shift
    ioff = 1
    call initialize_vshifte(ioff, nlon, nlat, wavetable, error_flag)

    ! Transfer regular grid values in (ureg, vreg) to offset grid in (uoff, voff)
    uoff = ZERO
    voff = ZERO

    write (stdout, 100) ioff, nlon, nlat

    call vshifte(ioff, nlon, nlat, uoff, voff, ureg, vreg, wavetable, error_flag)

    write (stdout, 200) error_flag

    if (error_flag == 0) then
        !     compute error in uoff, voff
        err2u = ZERO
        err2v = ZERO
        do j=1, nlon
            long = half_dlon+(j-1)*dlon
            sinp = sin(long)
            cosp = cos(long)
            do i=1, nlat
                lat = -HALF_PI + half_dlat + real(i - 1, kind=wp) * dlat
                sint = sin(lat)
                cost = cos(lat)
                x = cost*cosp
                y = cost*sinp
                z = sint
                ex = exp(x)
                ey = exp(y)
                ez = exp(z)
                emz = exp(-z)
                ue = -ex*sinp+emz*cost+ey*sint*sinp
                ve = ey*cosp+ez*cost-ex*sint*cosp
                err2u = err2u + (uoff(j, i)-ue)**2
                err2v = err2v + (voff(j, i)-ve)**2
            end do
        end do
        err2u = sqrt(err2u/(nlon*(nlat + 1)))
        err2v = sqrt(err2v/(nlon*(nlat + 1)))
        write (stdout, 300) err2u, err2v
    end if

    ! Release memory
    deallocate (wavetable)

end program testvshifte
