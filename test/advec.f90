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
!          d(phi)/dt = -(u,v) DOT gradient(phi)
!
!                    = -(u*gdphl + v*gdpht)
!
! ... required files
!
!     gradgc.f,shagc.f,shsgc.f,vhsgc.f,sphcom.f hrfft.f,gaqd.f
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
!     ncycle        exit number
!     time          model time in seconds
!     dt            time step
!     lambda        longitude
!     theta         latitude
!
!   the first dimension of the following two dimensional arrays
!   corresponds to the latitude index with values i=1,...,nlat
!   where i=1 is the northern most gaussian point thetag(i)
!   and i=nlat is the southern most gaussian point thetag(nlat).
!   the second dimension is longitude with values j=1,...,nlon
!   where j=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!
!     thetag(i)    vector of gaussian points on the full sphere which
!                  have north to south orientation as i=1,...,nlat
!
!     u(i,j)       east longitudinal velocity component
!     v(i,j)       latitudinal velocity component
!
!     phi(i,j)     the geopotential at t = time
!
!     phi_new(i,j)   the geopotential at t=time+dt
!
!     phi_old(i,j)   the geopotential at t=time-dt
!
!     gdphl(i,j)   the longitudinal derivative component of
!                  the gradient of phi
!
!                       gdphl = 1/(cos(theta))*d(phi)/dlambda

!
!     gdpht(i,j)   the latitudinal derivative component of
!                  the gradient of phi
!
!                       gdpht = d(phi)/dtheta

!
!   the following two dimensional arrays are nonzero in the triangle
!   n=1,...,nlat and m less than or equal to n.
!
!     ar(m,n),br(m,n)    spectral coefficients of phi
!
program advec
    implicit none
    !
    !     set grid size with parameter statements
    !
    integer nlat,nlon,lwork,ldwork,lvhsgc,lshagc

    !     parameter (nnlat=12,nnlon=23,ddt=1200.)
    integer, parameter :: nnlat=23
    integer, parameter :: nnlon=45
    real, parameter :: ddt=600.0
    !     parameter (nnlat=45,nnlon=90,ddt=300.)
    !
    !     set saved and unsaved work space lengths in terms of nnlat,nnlon
    !     (see documentation for shagc,vhsgc,vhsgci,gradgc for estimates)
    !
    integer, parameter :: nn15=nnlon+15
    integer, parameter :: llwork=4*nnlat*nnlon+2*nnlat*(nnlat+1)
    integer, parameter :: lldwork = 2*nnlat*(nnlat+1)+1
    integer, parameter :: llvhsgc = 7*nnlat*nnlat+nnlon+15
    integer, parameter :: llshagc = 5*nnlat*nnlat + nnlon+15
    !
    !     dimension arrays
    !
    real u(nnlat,nnlon),v(nnlat,nnlon)
    real phi_old(nnlat,nnlon),phi_new(nnlat,nnlon),phi(nnlat,nnlon)
    real exact_phi(nnlat,nnlon)
    real dpdt(nnlat,nnlon)
    real gdphl(nnlat,nnlon),gdpht(nnlat,nnlon), work(llwork)
    real dwork(lldwork)
    real wshagc(llshagc),wvhsgc(llvhsgc),wshsgc(llshagc)
    real ar(nnlat,nnlat),br(nnlat,nnlat)
    real thetag(nnlat),colat(nnlat)
    real dtheta(nnlat),dwts(nnlat)
    !
    !     set constants
    !
    real, parameter :: pi = acos(-1.0)
    real, parameter :: omega = (2.0*pi)/(12.0*24.0*3600.0)
    real, parameter :: p0 = 1000.0
    real, parameter :: re = 1.0/3
    real, parameter :: hzero = 1000.0
    real, parameter :: alphad = 60.0
    real, parameter :: alpha = pi*alphad/180
    real, parameter :: beta = pi/6
    !
    !     set one array and no equatorial symmetry
    !
    integer, parameter :: nt = 1
    integer, parameter :: isym = 0
    !
    !     set time step depending on resolution
    !
    real, parameter :: dt = ddt
    real, parameter :: tdt = dt+dt

    integer :: i, j, k, ierror, ier, error
    integer :: mprint, ncycle, ntime
    real :: errm, err2
    real :: p2, pmax, time, htime
    real :: xlm, sl, cl, st, ct, cthclh, cthslh, xlhat, clh, slh, cth, uhat, sth

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
    call gaqd(nlat,dtheta,dwts,dwork,ldwork,ier)

    thetag = 0.5*pi- dtheta
    colat = dtheta
    !
    !     preset saved work spaces for gradgc and shagc and shsgc
    !
    call vhsgci(nlat,nlon,wvhsgc,lvhsgc,dwork,ldwork,ierror)
    if (ierror /= 0) write(*,10) ierror
10  format(' error in vsgci = ',i5)
    call shagci(nlat,nlon,wshagc,lshagc,dwork,ldwork,ierror)
    if (ierror /= 0) write(*,20) ierror
20  format(' error in shagci = ',i5)
    call shsgci(nlat,nlon,wshsgc,lshagc,dwork,ldwork,ierror)
    if (ierror /= 0) write(*,21) ierror
21  format(' error in shsgci = ',i5)
    !
    !==> set vector velocities and cosine bell in geopotential
    !
    associate( &
        ca => cos(alpha), &
        sa => sin(alpha), &
        dlon => (pi+pi)/nlon &
        )
        do j=1,nlon
            xlm = (j-1)*dlon
            sl = sin(xlm)
            cl = cos(xlm)
            do i=1,nlat
                st = cos(colat(i))
                ct = sin(colat(i))
                sth = ca*st+sa*ct*cl
                cthclh = ca*ct*cl-sa*st
                cthslh = ct*sl
                xlhat = atanxy(cthclh,cthslh)
                clh = cos(xlhat)
                slh = sin(xlhat)
                cth = clh*cthclh+slh*cthslh
                uhat = omega*cth
                u(i,j) = (ca*sl*slh+cl*clh)*uhat
                v(i,j) = (ca*st*cl*slh-st*sl*clh+sa*ct*slh)*uhat
            end do
        end do
    end associate
    !
    !==> compute geopotential at t=-dt in phi_old and at t=0.0 in phi
    !    to start up leapfrog scheme
    !
    call gpot(-dt, alpha, beta, omega, hzero, re, nlat, nlon, nlat, colat, phi_old)
    call gpot(0.0, alpha, beta, omega, hzero, re, nlat, nlon, nlat, colat, phi)
    !
    !==> smooth geopotential at t=-dt and t=0. by synthesizing after analysis
    !
    call shagc(nlat,nlon,isym,nt,phi_old,nlat,nlon,ar,br,nlat, &
        nlat,wshagc,lshagc,work,lwork,ierror)
    if ( ierror /=0) write(*,26) ierror
    call shsgc(nlat,nlon,isym,nt,phi_old,nlat,nlon,ar,br,nlat, &
        nlat,wshsgc,lshagc,work,lwork,ierror)
    if ( ierror /=0) write(*,28) ierror
    call shagc(nlat,nlon,isym,nt,phi,nlat,nlon,ar,br,nlat, &
        nlat,wshagc,lshagc,work,lwork,ierror)
    if ( ierror /=0) write(*,26) ierror
    call shsgc(nlat,nlon,isym,nt,phi,nlat,nlon,ar,br,nlat, &
        nlat,wshsgc,lshagc,work,lwork,ierror)
    if ( ierror /=0) write(*,28) ierror
28  format(' ierror in shsgc = ',i5)
    !
    !==> compute l2 and max norms of geopotential at t=0.
    !
    p2 = norm2(phi)
    pmax = maxval(abs(phi))
    !
    !==> set number of time steps for 12 days
    !    (time to circumvent the earth)
    !
    ntime = int((12.0*24.0*3600.0)/dt+0.5)
    mprint = ntime/12
    time = 0.0
    ncycle = 0

    time_loop: do k=1,ntime+1
        !
        !       compute harmonic coefficients for phi at current time
        !
        call shagc(nlat,nlon,isym,nt,phi,nlat,nlon,ar,br,nlat, &
            nlat,wshagc,lshagc,work,lwork,ierror)
        if ( error /=0) write(*,26) ierror
26      format(' ierror in shagc = ',i5)
        !
        !       compute gradient of phi at current time
        !
        call gradgc(nlat,nlon,isym,nt,gdpht,gdphl,nlat,nlon,ar,br, &
            nlat,nlat,wvhsgc,lvhsgc,work,lwork,ierror)
        if ( error /=0) write(*,27) ierror
27      format(' ierror in gradgc = ',i5)
        !
        !==> compute the time derivative of phi, note that the sign
        !    of the last term is positive because the gradient is
        !    computed with respect to colatitude rather than latitude.
        !
        dpdt = -u*gdphl + v*gdpht

        if (mod(ncycle,mprint) == 0) then
            !
            !==> write variables
            !
            call gpot(time,alpha,beta,omega,hzero,re,nlat,nlon,nlat, &
                colat,exact_phi)

            errm = maxval(abs(exact_phi-phi))/pmax
            err2 = norm2(exact_phi-phi)/p2
            htime = time/3600.0

            write(*,390) ncycle,htime,dt,nlat,nlon,omega,hzero, &
                alphad,errm,err2
390         format(//' advecting cosine bell, test case 2',/ &
                ,' exit number              ',i10 &
                ,' model time in  hours      ',f10.2/ &
                ,' time step in seconds      ',f10.0 &
                ,' number of latitudes       ',i10/ &
                ,' number of longitudes      ',i10 &
                ,' rotation rate        ',1pe15.6/ &
                ,' mean height          ',1pe15.6 &
                ,' tilt angle                ',0pf10.2/ &
                ,' max geopot. error    ',1pe15.6 &
                ,' RMS geopot. error    ',1pe15.6)

        end if

        time = time + dt
        ncycle = ncycle+1
        !
        !==> update phi_old,phi for next time step
        !
        phi_new = phi_old + tdt*dpdt
        phi_old = phi
        phi = phi_new
        !
        !==> end of time loop
        !
    end do time_loop

contains

    subroutine gpot(t,alpha,beta,omega,hzero,re,nlat,nlon,idim, colat,h)
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
        implicit none
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real,    intent (in) :: t
        real,    intent (in) :: alpha
        real,    intent (in) :: beta
        real,    intent (in) :: omega
        real,    intent (in) :: hzero
        real,    intent (in) :: re
        integer, intent (in) :: nlat
        integer, intent (in) :: nlon
        integer, intent (in) :: idim
        real,    intent (in) :: colat(nlat)
        real,    intent (out) ::  h(idim,nlon)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer :: i, j
        real, parameter :: pi = acos(-1.0)
        real, parameter :: tpi = 2.0*pi
        real :: xc, yc, zc, x1, y1, z1
        real :: lambda, theta, st, ct, sth ,cthclh
        real :: cthslh,lhat,clh,slh,cth ,that, r
        !----------------------------------------------------------------------

        associate( lambdc => omega*t )
            call stoc(1.0, beta, lambdc, xc, yc, zc)
        end associate


        associate( &
            ca => cos(alpha), &
            sa => sin(alpha), &
            dlon => tpi/nlon &
            )
            do j=1,nlon
                lambda = (j-1)*dlon
                associate( &
                    cl => cos(lambda), &
                    sl => sin(lambda) &
                    )
                    do i=1,nlat
                        theta = colat(i)
                        st = cos(theta)
                        ct = sin(theta)
                        sth = ca*st+sa*ct*cl
                        cthclh = ca*ct*cl-sa*st
                        cthslh = ct*sl
                        lhat = atanxy(cthclh,cthslh)
                        clh = cos(lhat)
                        slh = sin(lhat)
                        cth = clh*cthclh+slh*cthslh
                        that = atanxy(sth,cth)
                        call stoc(1.0, that, lhat, x1, y1, z1)
                        associate( dist => norm2([x1-xc, y1-yc, z1-zc]) )
                            h(i,j) = 0.0
                            if (dist >= re) then
                                cycle
                            end if
                            r = 2.0*asin(dist/2)
                        end associate
                        if (r >= re) then
                            cycle
                        end if
                        h(i,j) = hzero*0.5*(cos(r*pi/re)+1.)
                    end do
                end associate
            end do
        end associate

    end subroutine gpot


    pure function atanxy(x,y) result (return_value)
        implicit none
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real, intent (in) :: x
        real, intent (in) :: y
        real              :: return_value
        !----------------------------------------------------------------------

        return_value = 0.0

        if (x == 0.0 .and. y == 0.0) then
            return
        end if

        return_value = atan2(y,x)

    end function atanxy


    pure subroutine ctos(x,y,z,r,theta,phi)
        implicit none
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real, intent (in) :: x, y, z
        real, intent (out) :: theta, r, phi
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real :: r1
        !----------------------------------------------------------------------

        r1 = norm2([x, y])

        if (r1 == 0.0) then
            phi = 0.
            theta = 0.
            if (z < 0.) then
                theta = acos(-1.0)
            end if
        else
            r = sqrt(r1+z**2)
            r1 = sqrt(r1)
            phi = atan2(y,x)
            theta = atan2(r1,z)
        end if

    end subroutine ctos


    pure subroutine stoc(r, theta, phi, x, y, z)
        implicit none
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real, intent (in) :: r
        real, intent (in) :: theta
        real, intent (in) :: phi
        real, intent (out) :: x, y, z
        !----------------------------------------------------------------------

        associate( st => sin(theta) )
            x = r*st*cos(phi)
            y = r*st*sin(phi)
            z = r*cos(theta)
        end associate

    end subroutine stoc

end program advec


