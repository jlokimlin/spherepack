!
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!  .                                                             .
!  .                  copyright (c) 1998 by UCAR                 .
!  .                                                             .
!  .       University Corporation for Atmospheric Research       .
!  .                                                             .
!  .                      all rights reserved                    .
!  .                                                             .
!  .                                                             .
!  .                          SPHEREPACK                         .
!  .                                                             .
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
!
!
! ... file advec.f
!
!     program advec solves the time-dependent linear advection 
!     equation for geopotential phi using the SPHEREPACK software
!
!          d(phi)/dt = -(u, v) DOT gradient(phi)
!
!                    = -(u*gdphl + v*gdpht)
!
! ... required files
!
!     gradgc.f, shagc.f, shsgc.f, vhsgc.f, sphcom.f hrfft.f, gaqd.f
!
!
! definitions:
!
!
!     nlat          number of gaussian latitudes excluding poles
!     nlon          number of distinct longitudes
!     omega         rotation rate of earth in radians per second
!     alpha         angle between axis of rotation and the coordinate
!                   axis
!     beta          latitude of the cosine bell
!     aa            radius of earth in meters
!     ncycle        cycle number
!     time          model time in seconds
!     dt            time step
!     lambda        longitude
!     theta         latitude
!
!   the first dimension of the following two dimensional arrays
!   corresponds to the latitude index with values i=1, ..., nlat
!   where i=1 is the northern most gaussian point thetag(i)
!   and i=nlat is the southern most gaussian point thetag(nlat).
!   the second dimension is longitude with values j=1, ..., nlon
!   where j=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!
!     thetag(i)    vector of gaussian points on the full sphere which
!                  have north to south orientation as i=1, ..., nlat
!
!     u(i, j)       east longitudinal velocity component
!     v(i, j)       latitudinal velocity component
!
!     phi(i, j)     the geopotential at t = time
!
!     phnew(i, j)   the geopotential at t=time+dt
!
!     phold(i, j)   the geopotential at t=time-dt
!
!     gdphl(i, j)   the longitudinal derivative component of
!                  the gradient of phi
!
!                       gdphl = 1/(cos(theta))*d(phi)/dlambda

!
!     gdpht(i, j)   the latitudinal derivative component of
!                  the gradient of phi
!
!                       gdpht = d(phi)/dtheta

!
!   the following two dimensional arrays are nonzero in the triangle
!   n=1, ..., nlat and m less than or equal to n.
!
!     ar(m, n), br(m, n)    spectral coefficients of phi
!
program advec
    use spherepack
    implicit none
    real(wp) :: alpha
    real(wp) :: alphad
    real(wp) :: beta
    real(wp) :: ca
    real(wp) :: cl
    real(wp) :: clh
    real(wp) :: ct
    real(wp) :: cth
    real(wp) :: cthclh
    real(wp) :: cthslh
    real(wp) :: ddt
    real(wp) :: dlon
    real(wp) :: dt
    real(wp) :: err2
    real(wp) :: errm
    real(wp) :: error
    real(wp) :: htime
    real(wp) :: hzero
    integer(ip) :: i
    integer(ip) :: ier
    integer(ip) :: ierror
    integer(ip) :: isym
    integer(ip) :: j
    integer(ip) :: k
    integer(ip) :: mprint
    integer(ip) :: ncycle
    integer(ip) :: nt
    integer(ip) :: ntime
    real(wp) :: omega
    real(wp) :: p0
    real(wp) :: p2

    real(wp) :: pmax
    real(wp) :: re
    real(wp) :: sa
    real(wp) :: sl
    real(wp) :: slh
    real(wp) :: st
    real(wp) :: sth
    real(wp) :: tdt
    real(wp) :: time
    real(wp) :: uhat
    real(wp) :: xlhat
    real(wp) :: xlm
    !
    !     set grid size with parameter statements
    !
    integer nnlat, nnlon, nn15, llwork, lldwork, llvhsgc, llshagc
    integer nlat, nlon, lwork, ldwork, lvhsgc, lshagc

    !     parameter (nnlat=12, nnlon=23, ddt=1200.)
    parameter (nnlat=23, nnlon=45, ddt=600.)
    !     parameter (nnlat=45, nnlon=90, ddt=300.)
    !
    !     set saved and unsaved work space lengths in terms of nnlat, nnlon
    !     (see documentation for shagc, vhsgc, vhsgci, gradgc for estimates)
    !
    parameter (nn15=nnlon+15)
    parameter(llwork=4*nnlat*nnlon+2*nnlat*(nnlat+1))
    parameter (lldwork = 2*nnlat*(nnlat+1)+1 )
    parameter (llvhsgc = 7*nnlat*nnlat+nnlon+15)
    parameter (llshagc = 5*nnlat*nnlat + nnlon+15)
    !
    !     dimension arrays
    !
    real u(nnlat, nnlon), v(nnlat, nnlon)
    real phold(nnlat, nnlon), phnew(nnlat, nnlon), phi(nnlat, nnlon)
    real pexact(nnlat, nnlon)
    real dpdt(nnlat, nnlon)
    real gdphl(nnlat, nnlon), gdpht(nnlat, nnlon), work(llwork)
    real dwork(lldwork)
    real wshagc(llshagc), wvhsgc(llvhsgc), wshsgc(llshagc)
    real ar(nnlat, nnlat), br(nnlat, nnlat)
    real thetag(nnlat), colat(nnlat)
    real dtheta(nnlat), dwts(nnlat)
    !
    !     set constants
    !
    omega = (pi+pi)/(12.*24.*3600.)
    p0 = 1000.
    re = 1.0/3.0
    hzero = 1000.
    alphad = 60.
    alpha = pi*alphad/180.
    beta = pi/6.
    !
    !     set one array and no equatorial symmetry
    !
    nt = 1
    isym = 0
    !
    !     set time step depending on resolution
    !
    dt = ddt
    tdt = dt+dt
    !
    !     set work space length arguments
    !
    lwork = llwork
    ldwork = lldwork
    lshagc = llshagc
    lvhsgc = llvhsgc
    !
    !     set grid size arguments
    !
    nlat = nnlat
    nlon = nnlon
    !
    !     compute nlat latitudinal gaussian points in thetag  with
    !     north to south orientation using gaqd from SPHEREPACK
    !
    call compute_gaussian_latitudes_and_weights(nlat, dtheta, dwts, ier)
    do  i=1, nlat
        thetag(i) = 0.5*pi- dtheta(i)
        colat(i) = dtheta(i)
    end do
    !
    !     preset saved work spaces for gradgc and shagc and shsgc
    !
    call vhsgci(nlat, nlon, wvhsgc, lvhsgc, dwork, ldwork, ierror)
    if(ierror /= 0) write(*, 10) ierror
10  format(' error in vsgci = ', i5)
    call shagci(nlat, nlon, wshagc, ierror)
    if(ierror /= 0) write(*, 20) ierror
20  format(' error in shagci = ', i5)
    call shsgci(nlat, nlon, wshsgc, lshagc, dwork, ldwork, ierror)
    if(ierror /= 0) write(*, 21) ierror
21  format(' error in shsgci = ', i5)
    !
    !     set vector velocities and cosine bell in geopotential
    !
    ca = cos(alpha)
    sa = sin(alpha)
    dlon = (pi+pi)/nlon
    do j=1, nlon
        xlm = (j-1)*dlon
        sl = sin(xlm)
        cl = cos(xlm)
        do i=1, nlat
            st = cos(colat(i))
            ct = sin(colat(i))
            sth = ca*st+sa*ct*cl
            cthclh = ca*ct*cl-sa*st
            cthslh = ct*sl
            xlhat = atanxy(cthclh, cthslh)
            clh = cos(xlhat)
            slh = sin(xlhat)
            cth = clh*cthclh+slh*cthslh
            uhat = omega*cth
            u(i, j) = (ca*sl*slh+cl*clh)*uhat
            v(i, j) = (ca*st*cl*slh-st*sl*clh+sa*ct*slh)*uhat
        end do
    end do
    !
    !       compute geopotential at t=-dt in phold and at t=0.0 in phi
    !       to start up leapfrog scheme
    !
    call gpot(-dt, alpha, beta, omega, hzero, re, nlat, nlon, &
        nlat, colat, phold)
    call gpot(0., alpha, beta, omega, hzero, re, nlat, nlon, &
        nlat, colat, phi)
    !
    !     smooth geopotential at t=-dt and t=0. by synthesizing after analysis
    !
    call shagc(nlat, nlon, isym, nt, phold, nlat, nlon, ar, br, nlat, &
        nlat, wshagc, lshagc, work, lwork, ierror)
    if ( ierror /=0) write(*, 26) ierror
    call shsgc(nlat, nlon, isym, nt, phold, nlat, nlon, ar, br, nlat, &
        nlat, wshsgc, lshagc, work, lwork, ierror)
    if ( ierror /=0) write(*, 28) ierror
    call shagc(nlat, nlon, isym, nt, phi, nlat, nlon, ar, br, nlat, &
        nlat, wshagc, lshagc, work, lwork, ierror)
    if ( ierror /=0) write(*, 26) ierror
    call shsgc(nlat, nlon, isym, nt, phi, nlat, nlon, ar, br, nlat, &
        nlat, wshsgc, lshagc, work, lwork, ierror)
    if ( ierror /=0) write(*, 28) ierror
28  format(' ierror in shsgc = ', i5)
    !
    !     compute l2 and max norms of geopotential at t=0.
    !
    p2 = 0.0
    pmax = 0.0
    do j=1, nlon
        do i=1, nlat
            pmax = amax1(abs(phi(i, j)), pmax)
            p2 = p2 + phi(i, j)**2
        end do
    end do
    p2 = sqrt(p2)
    !
    !     set number of time steps for 12 days
    !     (time to circumvent the earth)
    !
    ntime = int((12.*24.*3600.)/dt+0.5)
    mprint = ntime/12
    time = 0.0
    ncycle = 0
    do k=1, ntime+1
        !
        !       compute harmonic coefficients for phi at current time
        !
        call shagc(nlat, nlon, isym, nt, phi, nlat, nlon, ar, br, nlat, &
            nlat, wshagc, lshagc, work, lwork, ierror)
        if ( error /=0) write(*, 26) ierror
26      format(' ierror in shagc = ', i5)
        !
        !       compute gradient of phi at current time
        !
        call gradgc(nlat, nlon, isym, nt, gdpht, gdphl, nlat, nlon, ar, br, &
            nlat, nlat, wvhsgc, lvhsgc, work, lwork, ierror)
        if ( error /=0) write(*, 27) ierror
27      format(' ierror in gradgc = ', i5)
                !
                !       compute the time derivative of phi, note that the sign
                !       of the last term is positive because the gradient is
                !       computed with respect to colatitude rather than latitude.
                !
        do j=1, nlon
            do i=1, nlat
                dpdt(i, j) = -u(i, j)*gdphl(i, j) + v(i, j)*gdpht(i, j)
            end do
        end do
        !
        if (mod(ncycle, mprint) == 0) then
            !
            !     write variables
            !
            err2 = 0.0
            errm = 0.0
            call gpot(time, alpha, beta, omega, hzero, re, nlat, nlon, nlat, &
                colat, pexact)
            do j=1, nlon
                do i=1, nlat
                    err2 = err2 + (pexact(i, j)-phi(i, j))**2
                    errm = amax1(abs(pexact(i, j)-phi(i, j)), errm)
                end do
            end do
            errm = errm/pmax
            err2 = sqrt(err2)/p2
            htime = time/3600.
            write(*, 390) ncycle, htime, dt, nlat, nlon, omega, hzero, &
                alphad, errm, err2
390         format(//' advecting cosine bell, test case 2', / &
                , ' cycle number              ', i10 &
                , ' model time in  hours      ', f10.2/ &
                , ' time step in seconds      ', f10.0 &
                , ' number of latitudes       ', i10/ &
                , ' number of longitudes      ', i10 &
                , ' rotation rate        ', 1pe15.6/ &
                , ' mean height          ', 1pe15.6 &
                , ' tilt angle                ', 0pf10.2/ &
                , ' max geopot. error    ', 1pe15.6 &
                , ' RMS geopot. error    ', 1pe15.6)

        end if
        time = time + dt
        ncycle = ncycle+1
            !
            !       update phold, phi for next time step
            !
        do j=1, nlon
            do i=1, nlat
                phnew(i, j) = phold(i, j) + tdt*dpdt(i, j)
                phold(i, j) = phi(i, j)
                phi(i, j) = phnew(i, j)
            end do
        end do
    !
    !     end of time loop
    !
    end do

contains

    subroutine gpot(t, alpha, beta, omega, hzero, re, nlat, nlon, idim, &
        colat, h)
        implicit none
        real(wp) :: alpha
        real(wp) :: beta
        real(wp) :: ca
        real(wp) :: cl
        real(wp) :: clh
        real(wp) :: colat
        real(wp) :: ct
        real(wp) :: cth
        real(wp) :: cthclh
        real(wp) :: cthslh
        real(wp) :: dist
        real(wp) :: dlon
        real(wp) :: h
        real(wp) :: hzero
        integer(ip) :: i
        integer(ip) :: idim
        integer(ip) :: j
        integer(ip) :: nlat
        integer(ip) :: nlon
        real(wp) :: omega
        real(wp) :: pi
        real(wp) :: r
        real(wp) :: re
        real(wp) :: sa
        real(wp) :: sl
        real(wp) :: slh
        real(wp) :: st
        real(wp) :: sth
        real(wp) :: t
        real(wp) :: that
        real(wp) :: theta
        real(wp) :: tpi
        real(wp) :: x1
        real(wp) :: xc
        real(wp) :: y1
        real(wp) :: yc
        real(wp) :: z1
        real(wp) :: zc
        !
        !     computes advecting cosine bell on a tilted grid a time t.
        !
        ! input parameters
        !
        !     t      time in seconds
        !
        !     alpha  tilt angle in radians
        !
        !     beta   colatitude of cosine bell in untilted coordinate
        !            system in radians
        !
        !     omega  angular velocity in radians per second
        !
        !     hzero  maximum value of cosine bell
        !
        !     re     radius of support for cosine bell in radians
        !
        !     nlat   number of latitudes including the poles
        !
        !     nlon   number of distinct longitude lines
        !
        !     idim   first dimension of output array h
        !
        !     colat  vector of Gauss colatitude grid points
        !
        ! output parameter
        !
        !     h      an nlat by nlon array containing the geopotential
        !
        !             on a tilted grid
        !
        dimension h(idim, nlon), colat(nlat)
        real lambda, lambdc, lhat
        lambdc = omega*t
        call stoc(1., beta, lambdc, xc, yc, zc)
        ca = cos(alpha)
        sa = sin(alpha)
        pi = 4.*atan(1.)
        tpi = pi+pi
        dlon = tpi/nlon
        do 10 j=1, nlon
            lambda = (j-1)*dlon
            cl = cos(lambda)
            sl = sin(lambda)
            do 10 i=1, nlat
                theta = colat(i)
                st = cos(theta)
                ct = sin(theta)
                sth = ca*st+sa*ct*cl
                cthclh = ca*ct*cl-sa*st
                cthslh = ct*sl
                lhat = atanxy(cthclh, cthslh)
                clh = cos(lhat)
                slh = sin(lhat)
                cth = clh*cthclh+slh*cthslh
                that = atanxy(sth, cth)
                call stoc(1., that, lhat, x1, y1, z1)
                dist = sqrt((x1-xc)**2+(y1-yc)**2 &
                    +(z1-zc)**2)
                h(i, j) = 0.
                if(dist >= re) go to 10
                r = 2.*asin(dist/2.)
                if(r >= re) go to 10
                h(i, j) = hzero*.5*(cos(r*pi/re)+1.)
10          continue

        end subroutine gpot


        real function atanxy(x, y)
            implicit none
            real(wp) :: x
            real(wp) :: y
            atanxy = 0.
            if(x==0. .and. y==0.) return
            atanxy = atan2(y, x)

        end function atanxy


        subroutine ctos(x, y, z, r, theta, phi)
            implicit none
            real(wp) :: phi
            real(wp) :: r
            real(wp) :: r1
            real(wp) :: theta
            real(wp) :: x
            real(wp) :: y
            real(wp) :: z
            r1 = x*x+y*y
            if(r1 /= 0.) go to 10
            phi = 0.
            theta = 0.
            if(z < 0.) theta = 4.*atan(1.)
            return
10          r = sqrt(r1+z*z)
            r1 = sqrt(r1)
            phi = atan2(y, x)
            theta = atan2(r1, z)

        end subroutine ctos


        subroutine stoc(r, theta, phi, x, y, z)
            implicit none
            real(wp) :: phi
            real(wp) :: r
            real(wp) :: st
            real(wp) :: theta
            real(wp) :: x
            real(wp) :: y
            real(wp) :: z
            st = sin(theta)
            x = r*st*cos(phi)
            y = r*st*sin(phi)
            z = r*cos(theta)

        end subroutine stoc

    end program advec
