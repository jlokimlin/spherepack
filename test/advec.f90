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

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        GaussianSphere

    ! Explicit typing only
    implicit none

    !----------------------------------------------------------------------
    ! Dictionary: global variables
    !----------------------------------------------------------------------
    real (wp),    parameter :: PI = acos(-1.0_wp)
    real (wp),    parameter :: DEGREES_TO_RADIANS = PI/180
    integer (ip), parameter :: TIME_TO_CIRCUMVENT_THE_EARTH = 12*24*3600
    !----------------------------------------------------------------------


    ! Declare derived data type
    type, extends (GaussianSphere) :: Solver
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        real (wp) :: TIME_STEP = 6.0e+2_wp
        real (wp) :: ROTATION_RATE_OF_EARTH = 2.0_wp * PI/TIME_TO_CIRCUMVENT_THE_EARTH
        real (wp) :: LATITUDE_OF_COSINE_BELL = PI/6
        real (wp) :: RADIUS_OF_EARTH_IN_METERS = 1.0_wp/3
        real (wp) :: TILT_ANGLE_IN_DEGREES = 60.0_wp
        real (wp) :: TILT_ANGLE = 60.0_wp * DEGREES_TO_RADIANS
        real (wp) :: MAXIMUM_VALUE_OF_COSINE_BELL = 1.0e+3_wp
        !----------------------------------------------------------------------
    end type Solver



    call advecting_cosine_bell()



contains



    subroutine advecting_cosine_bell()
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        type (Solver)                 :: solver_dat
        integer (ip), parameter       :: NLONS = 45
        integer (ip), parameter       :: NLATS = 23
        integer (ip)                  :: i, j !! Counters
        integer (ip)                  :: mprint, ncycle, ntime
        real (wp)                     :: u(NLATS,NLONS),v(NLATS,NLONS)
        real (wp)                     :: phi_old(NLATS,NLONS),phi_new(NLATS,NLONS)
        real (wp)                     :: phi(NLATS,NLONS)
        real (wp)                     :: exact_phi(NLATS,NLONS)
        real (wp)                     :: dpdt(NLATS,NLONS)
        real (wp)                     :: gdphl(NLATS,NLONS),gdpht(NLATS,NLONS)
        real (wp)                     :: p0_l2, p0_max, time, htime
        real (wp)                     :: xlhat, uhat, cth
        character (len=*), parameter  :: write_format = &
            '( A,i10 ,A,f10.2/ ,A,f10.0 ,A,i10/ ,A,i10 '&
            //',A,1pe15.6/,A,1pe15.6 ,A,0pf10.2/ A,1pe15.6 ,A,1pe15.6)'
        !----------------------------------------------------------------------

        !
        !==> Set up workspace arrays
        !
        call solver_dat%create(nlat=NLATS, nlon=NLONS)
        !
        !==> set vector velocities and cosine bell in geopotential
        !
        associate( &
            colat => solver_dat%grid%latitudes, &
            omega => solver_dat%ROTATION_RATE_OF_EARTH, &
            alpha => solver_dat%TILT_ANGLE &
            )
            associate( &
                cosa => cos(alpha), &
                sina => sin(alpha) &
                )
                do j=1,NLONS
                    associate( xlm => solver_dat%grid%longitudes(j) )
                        associate( &
                            sinp => sin(xlm), &
                            cosp => cos(xlm) &
                            )
                            do i=1,NLATS
                                associate( &
                                    cost => cos(colat(i)), &
                                    sint => sin(colat(i)) &
                                    )
                                    associate( &
                                        sth => cosa*cost+sina*sint*cosp, &
                                        cthclh => cosa*sint*cosp-sina*cost, &
                                        cthslh => sint*sinp &
                                        )
                                        xlhat = atanxy(cthclh,cthslh)
                                        associate( &
                                            cosl => cos(xlhat), &
                                            sinl => sin(xlhat) &
                                            )
                                            cth = cosl*cthclh+sinl*cthslh
                                            uhat = omega*cth
                                            u(i,j) = (cosa*sinp*sinl+cosp*cosl)*uhat
                                            v(i,j) = &
                                                (cosa*cost*cosp*sinl &
                                                -cost*sinp*cosl+sina*sint*sinl)*uhat
                                        end associate
                                    end associate
                                end associate
                            end do
                        end associate
                    end associate
                end do
            end associate
        end associate
        !
        !==> compute geopotential at t=-dt in phi_old and at t=0.0 in phi
        !    to start up leapfrog scheme
        !
        associate( dt => solver_dat%TIME_STEP )
            call get_geopotential(solver_dat, -dt, phi_old)
            call get_geopotential(solver_dat, 0.0_wp, phi)
        end associate
        !
        !==> smooth geopotential at t=-dt and t=0. by synthesizing after analysis
        !
        call solver_dat%perform_scalar_analysis(phi_old)
        call solver_dat%perform_scalar_synthesis(phi_old)
        call solver_dat%perform_scalar_analysis(phi)
        call solver_dat%perform_scalar_synthesis(phi)
        !
        !==> compute l2 and max norms of geopotential at t=0.
        !
        p0_l2 = norm2(phi)
        p0_max = maxval(abs(phi))
        !
        !==> set number of time steps for 12 days
        !    (time to circumvent the earth)
        !
        associate( dt => solver_dat%TIME_STEP )
            ntime = int(real(TIME_TO_CIRCUMVENT_THE_EARTH, kind=wp)/dt + 0.5_wp, kind=ip)
            mprint = ntime/12
            time = 0.0_wp
            ncycle = 0
        end associate
        !
        !==> Start time loop
        !
        time_loop: do while (ncycle <= ntime)
            !
            !==> Compute gradient of phi at current time
            !
            call solver_dat%get_gradient(phi, gdpht, gdphl)
            !
            !==> Compute the time derivative of phi, note that the sign
            !    of the last term is positive because the gradient is
            !    computed with respect to colatitude rather than latitude.
            !
            dpdt = -u * gdphl + v * gdpht
            !
            !==> write variables to standard output
            !
            if (mod(ncycle, mprint) == 0) then
                !
                !==> Compute exact solution
                !
                call get_geopotential(solver_dat, time, exact_phi)
                !
                !==> Compute errors
                !
                associate( &
                    htime => time/3600, &
                    dt => solver_dat%TIME_STEP, &
                    nlats => solver_dat%NUMBER_OF_LATITUDES, &
                    nlons => solver_dat%NUMBER_OF_LONGITUDES, &
                    omega => solver_dat%ROTATION_RATE_OF_EARTH, &
                    hzero => solver_dat%MAXIMUM_VALUE_OF_COSINE_BELL, &
                    alphad => solver_dat%TILT_ANGLE_IN_DEGREES, &
                    errm => maxval(abs(exact_phi-phi))/p0_max, &
                    err2 => norm2(exact_phi-phi)/p0_l2 &
                    )
                    write( stdout, '(A)') ''
                    write( stdout, '(A)' ) ' advecting cosine bell, test case 2'
                    write( stdout, fmt = write_format ) &
                        ' exit number              ', ncycle, &
                        '  model time in  hours      ', htime, &
                        ' time step in seconds      ', dt, &
                        ' number of latitudes       ', nlats, &
                        ' number of longitudes      ', nlons, &
                        ' rotation rate        ', omega, &
                        ' mean height          ', hzero, &
                        ' tilt angle                ', alphad, &
                        ' max geopot. error    ', errm, &
                        ' RMS geopot. error    ', err2
                end associate
            end if

            !
            !==> Update various quantities
            !
            associate( dt => solver_dat%TIME_STEP )
                !
                !==> increment
                !
                time = time + dt
                ncycle = ncycle + 1
                !
                !==> update phi_old,phi for next time step
                !
                phi_new = phi_old + 2.0_wp * dt * dpdt
                phi_old = phi
                phi = phi_new
            end associate
            !
            !==> end of time loop
            !
        end do time_loop

    end subroutine advecting_cosine_bell



    subroutine get_geopotential(this, t, geopot)
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
        !     geopot     an nlat by nlon array containing the geopotential
        !
        !             on a tilted grid
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        class (Solver), intent (in out) :: this
        real (wp),      intent (in)     :: t
        real (wp),      intent (out)    :: geopot(:,:)
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: i, j
        real (wp)    :: xc, yc, zc, x1, y1, z1
        real (wp)    :: lambda, theta, cost, sint, sth ,cthclh
        real (wp)    :: cthslh,lhat,cosl,sinl,cth ,that, r
        !----------------------------------------------------------------------

        associate( &
            nlat => this%NUMBER_OF_LATITUDES, &
            nlon => this%NUMBER_OF_LONGITUDES, &
            omega => this%ROTATION_RATE_OF_EARTH, &
            beta => this%LATITUDE_OF_COSINE_BELL, &
            re => this%RADIUS_OF_EARTH_IN_METERS, &
            alpha => this%TILT_ANGLE, &
            hzero => this%MAXIMUM_VALUE_OF_COSINE_BELL,&
            colat => this%grid%latitudes &
            )

            associate( lambdc => omega*t )
                call sph2cart(1.0, beta, lambdc, xc, yc, zc)
            end associate


            associate( &
                cosa => cos(alpha), &
                sina => sin(alpha) &
                )
                do j=1,nlon
                    lambda = this%grid%longitudes(j)
                    associate( &
                        cosp => cos(lambda), &
                        sinp => sin(lambda) &
                        )
                        do i=1,nlat
                            theta = colat(i)
                            cost = cos(theta)
                            sint = sin(theta)
                            sth = cosa*cost+sina*sint*cosp
                            cthclh = cosa*sint*cosp-sina*cost
                            cthslh = sint*sinp
                            lhat = atanxy(cthclh,cthslh)
                            cosl = cos(lhat)
                            sinl = sin(lhat)
                            cth = cosl*cthclh+sinl*cthslh
                            that = atanxy(sth,cth)
                            !
                            !== Compute scaled radial unit vector
                            !
                            call sph2cart(1.0_wp, that, lhat, x1, y1, z1)
                            !
                            !==> compute distance
                            !
                            associate( dist => norm2([x1-xc, y1-yc, z1-zc]) )
                                !
                                !==> Initialize geopotential
                                !
                                geopot(i,j) = 0.0_wp
                                if (dist >= re) then
                                    cycle
                                end if
                                r = 2.0_wp * asin(dist/2)
                            end associate
                            if (r >= re) then
                                cycle
                            end if
                            !
                            !==> Set geopotential
                            geopot(i,j) = hzero * (cos(r*PI/re)+1.0_wp)/2
                        end do
                    end associate
                end do
            end associate
        end associate

    end subroutine get_geopotential



    pure function atanxy(x,y) result (return_value)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp), intent (in) :: x
        real (wp), intent (in) :: y
        real                   :: return_value
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (wp), parameter :: ZERO = nearest(1.0_wp, 1.0_wp)-nearest(1.0_wp, -1.0_wp)
        !--------------------------------------------------------------------------------

        return_value = 0.0_wp

        if (x == ZERO .and. y == ZERO) then
            return
        end if

        return_value = atan2(y,x)

    end function atanxy



    pure subroutine sph2cart(radius, theta, phi, x, y, z)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        real (wp), intent (in)  :: radius
        real (wp), intent (in)  :: theta
        real (wp), intent (in)  :: phi
        real (wp), intent (out) :: x, y, z
        !----------------------------------------------------------------------

        x = radius * sin(theta) * cos(phi)
        y = radius * sin(theta) * sin(phi)
        z = radius * cos(theta)

    end subroutine sph2cart



end program advec
